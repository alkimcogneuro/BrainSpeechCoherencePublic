function [SpeechEEGCoherence_Results] = SpeechEEGCoherence_ShuffleTest(eeg_struct, speech_rawdata, band_low, band_high, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % SHUFFLE-BASED SIGNIFICANCE TEST FOR BAND-AVERAGED SPEECH-EEG MAGNITUDE SQUARED COHERENCE
    % ===============================================================================================================
    % Tests, for one subject, at every EEG channel, whether the TRUE (correctly-paired) trial-averaged,
    % BAND-AVERAGED magnitude squared coherence (MSC) between the speech envelope and the EEG is
    % greater than expected by chance prepped
    % We perform a separate permutation test at each channel, 
    % run on the same quantity that BandAverageMSC.m produces, rather than one test per channel
    % per frequency bin (CSD and PSD are still computed at full frequency resolution, but the final test is on the band-averaged MSC).
    % 
    % "Expected by chance" is defined as follows: 
    % 1. Take a random derangement of the trial order (derange.m -- a permutation with
    %       no trial mapped to itself), so every EEG trial is analyzed against the speech epoch belonging to a different trial,
    % 2. run BandAverageMSC.m on that shuffled pairing, so you get a single band-averaged MSC value per channel for that shuffle.
    % 3. Repeat 1-2 many times.  This builds a null distribution of band-averaged MSC values, per channel, 
    %       against which the true band-averaged MSC is compared.
    %
    % This function calls PreprocessSpeechEEGForCoherence.m once to do the expensive, pairing-independent preprocessing
    %  (z-scoring, band-pass filtering every trial/channel, and extracting & filtering every trial's speech
    %   envelope)
    % This step also validates/determines the Welch analysis-window length 
    % used for every trial's CSD/PSD estimate, see options.WindowDurationSec below), 
    % 
    % After preprocessing, we call ComputeCoherenceForPairing.m once for the true pairing and once per shuffle, 
    % band-averaging each result. 
    % For full rationale, see PreprocessSpeechEEGForCoherence.m, ComputeCoherenceForPairing.m, and WelchCrossSpectralDensity.m.
    % 
    % For each trial, CSD/PSD come from WelchCrossSpectralDensity.m, which does
    % 1. calculates CSD in several overlapping windows within the trial, averaged together, rather than one window spanning the whole trial
    % 2. then averaged again across trials; still computed at full [Welch-window] frequency resolution
    % 
    % Averaging MSC values within the band of interest (band_low <= f <= band_high) is done later, by BandAverageMSC.m, 
    % which takes the full CSD/PSD spectra as input and returns a single band-averaged MSC value per channel.
    % 
    % Significance Test is one-sided: p(channel) = P(null band-MSC >= true band-MSC), using the standard
    % permutation-test estimator that never returns exactly zero:
    %     p = (1 + #{shuffles with null band-MSC >= true band-MSC}) / (NumShuffles + 1)
    % The smallest possible p-value is therefore 1/(NumShuffles+1), so
    % NumShuffles should be chosen with that resolution in mind (e.g. for 1000 shuffles, smallest
    % possible p is ~0.001).
    %
    % NOTE on multiple comparisons: this test is still run independently at every channel (e.g. 64
    % tests for a 64-channel montage), so we might consider FDR correction (e.g. Benjamini-Hochberg) across
    % channels before declaring significance, even though the frequency dimension has been collapsed.
    % At the moment, we are not using this test for any strong inference, 
    % so we are not correcting for multiple comparisons across channels, but that could be added later if desired.
    % 
    %   ARGUMENTS:
    %   eeg_struct, speech_rawdata : same meaning as in AnalyzeSpeechEEGCoherence_Dataset.m
    %   band_low, band_high : the frequency band (Hz) to average MSC over, passed straight through
    %     to BandAverageMSC.m (band_low <= band_high, inclusive on both ends).
    %   highpass_cutoff, lowpass_cutoff (defaults 2, 35) : same meaning as in
    %     AnalyzeSpeechEEGCoherence_Dataset.m -- the broadband filter applied before CSD estimation.
    %     (band_low/band_high select which part of that already-analyzed spectrum gets averaged;
    %     they are unrelated to highpass_cutoff/lowpass_cutoff, which control what goes into the FFT
    %     in the first place. band_low/band_high should fall within [highpass_cutoff, lowpass_cutoff].)
    %   options.StartTimeOffset, options.EpochDuration : same meaning as in
    %     AnalyzeSpeechEEGCoherence_Dataset.m
    %   options.NumShuffles (default 0) : number of random derangements to draw for the null
    %     distribution.
    %   options.RandomSeed (default NaN = don't seed) : if provided, calls rng(RandomSeed) before
    %     drawing shuffles, for reproducibility. NOTE: this reseeds MATLAB's global random stream.
    %   options.AvgMethod (default 'power-weighted') : passed through to BandAverageMSC.m -- either
    %     'power-weighted' (recommended, matches how CSD/PSD are already averaged across trials in
    %     this pipeline) or 'simple'. See BandAverageMSC.m for the rationale.
    %   options.WindowDurationSec, options.OverlapFraction, options.MinFreqHz,
    %     options.MinCyclesPerWindow : passed straight through to PreprocessSpeechEEGForCoherence.m
    %     to control the Welch analysis window used for every trial's CSD/PSD estimate -- see that
    %     file's header (and RecommendWelchWindowLength.m) for the full rationale and defaults.
    %     WindowDurationSec should be chosen deliberately per dataset once you know its trial
    %     length (e.g. via RecommendWelchWindowLength.m) rather than left to the auto floor for
    %     every dataset, since the right window length for a 30 sec trial is not the right window
    %     length for a much shorter one.
    %
    % Returns SpeechEEGCoherence_Results, a struct with fields:
    %   .True_MSC    [Num_channels x 1] true (observed) band-averaged MSC, per channel -- identical
    %                to what BandAverageMSC(CoherenceResults, band_low, band_high) would return for
    %                this same data, where CoherenceResults comes from
    %                AnalyzeSpeechEEGCoherence_Dataset with RandomizeOnsets=false.
    %   .PValues     [Num_channels x 1] one-sided permutation p-values, per channel.
    %   .ZScores     [Num_channels x 1] (True_MSC - Null_Mean) ./ Null_Std, per channel (a
    %                parametric-flavored companion to the p-values -- useful for ranking/plotting,
    %                but the p-values are the actual significance test).
    %   .Null_Mean, .Null_Std   [Num_channels x 1] mean/std of the null band-averaged MSC distribution.
    %   .Null_MSC    [NumShuffles x Num_channels] the full null distribution (band-averaged MSC is
    %                a single number per channel per shuffle, so this is small -- a few hundred KB
    %                even for thousands of shuffles -- and is always returned).
    %   .PairingsUsed [NumShuffles x Num_trials] the derangement used for each shuffle.
    %   .band_low, .band_high, .band_freqs   the requested band and the actual frequency bins
    %                (Hz) that fell inside it and were averaged (see BandAverageMSC.m).
    %   .AvgMethod, .NumShuffles, .nfft, .eeg_Fs, .highpass_cutoff, .lowpass_cutoff, .RandomSeed
    %   .label, .Chanlocs, .Condition, .Subj_id, .audio_file   (passthrough metadata)
    %   .TestDescription  (human-readable string describing the test, for provenance)
    % ===============================================================================================================
    arguments
        eeg_struct
        speech_rawdata
        band_low (1,1) double {mustBeNonnegative}  % the frequency band (Hz) to average MSC over, passed straight through to BandAverageMSC.m (band_low <= band_high, inclusive on both ends).
        band_high (1,1) double {mustBeNonnegative}  
        highpass_cutoff (1,1) double {mustBeReal} = 2
        lowpass_cutoff  (1,1) double {mustBeReal} = 35
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN
        options.EpochDuration (1,1) double {mustBeReal} = NaN
        options.NumShuffles (1,1) double {mustBeInteger} = 0
        options.RandomSeed (1,1) double {mustBeReal} = NaN
        options.AvgMethod (1,:) char = 'power-weighted'
        options.WindowDurationSec (1,1) double {mustBeReal} = NaN
        options.OverlapFraction (1,1) double {mustBeReal} = 0.5
        options.MinFreqHz (1,1) double {mustBeReal} = highpass_cutoff
        options.MinCyclesPerWindow (1,1) double {mustBeReal} = 3
    end

    % Make sure band_low < band_high, otherwise input to BandAverageMSC.m makes no sense.
    if band_low >= band_high
        error('ShuffleTest_SpeechEEGCoherence:InvalidBand', ...
            'band_low (%.2f Hz) must be less than band_high (%.2f Hz).', band_low, band_high);
    end

    fprintf('=== ShuffleTest_SpeechEEGCoherence: subject %s, band [%.2f %.2f] Hz, %d shuffles ===\n', ...
        eeg_struct.Subj_id, band_low, band_high, options.NumShuffles);

    if ~isnan(options.RandomSeed)
        fprintf(' -- Seeding global random stream with RandomSeed = %d for reproducibility.\n', options.RandomSeed);
        rng(options.RandomSeed);
    end

    % ------------------------------------------------------------------------------------------------------------
    % Step 1: Expensive, pairing-independent preprocessing, done just once for the dataset. 
    % Computed at full frequency resolution inside ComputeCoherenceForPairing below. 
    % Only the accumulation/statistics, after initial CSD/PSD/MSC computation, are band-averaged.
    % ------------------------------------------------------------------------------------------------------------
    EEG_Speech_prepped = PreprocessSpeechEEGForCoherence(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, ...
        StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration, ...
        WindowDurationSec=options.WindowDurationSec, OverlapFraction=options.OverlapFraction, ...
        MinFreqHz=options.MinFreqHz, MinCyclesPerWindow=options.MinCyclesPerWindow);

    Num_channels = EEG_Speech_prepped.Num_channels;
    Num_trials = EEG_Speech_prepped.Num_trials;

    % ------------------------------------------------------------------------------------------------------------
    % Step 2: compute true coherence for the observed pairing 
    % Each EEG trial is paired with its own speech,
    % we compute trial-averaged CSD/PSD for that pairing (ComputeCoherenceForPairing.m),
    % then band-average (BandAverageMSC.m) to get a single band-averaged MSC value per channel.
    % ------------------------------------------------------------------------------------------------------------
    fprintf(' -- Computing true (observed) band-averaged coherence...\n');
    TrueResults = ComputeCoherenceForPairing(EEG_Speech_prepped, 1:Num_trials);
    TrueBandResults = BandAverageMSC(TrueResults, band_low, band_high, AvgMethod=options.AvgMethod);
    True_MSC = TrueBandResults.msc_band;
    band_freqs = TrueBandResults.band_freqs;

    % Package some of the metadata and results into the output struct 
    % for return even if we skip the shuffle test below.
    % If we do run the shuffle test, then we information about the null distribution and p-values will be added to this struct.
    SpeechEEGCoherence_Results = struct();
    SpeechEEGCoherence_Results.True_MSC = True_MSC;
    SpeechEEGCoherence_Results.band_low = band_low;
    SpeechEEGCoherence_Results.band_high = band_high;
    SpeechEEGCoherence_Results.band_freqs = band_freqs;
    SpeechEEGCoherence_Results.AvgMethod = options.AvgMethod;
    SpeechEEGCoherence_Results.nfft = EEG_Speech_prepped.nfft;
    SpeechEEGCoherence_Results.WindowDurationSec = EEG_Speech_prepped.WelchWindowSamples / EEG_Speech_prepped.eeg_Fs;
    SpeechEEGCoherence_Results.OverlapFraction = EEG_Speech_prepped.WelchOverlapFraction;
    SpeechEEGCoherence_Results.NumSegmentsPerTrial = EEG_Speech_prepped.WelchNumSegmentsPerTrial;
    SpeechEEGCoherence_Results.eeg_Fs = EEG_Speech_prepped.eeg_Fs;
    SpeechEEGCoherence_Results.highpass_cutoff = highpass_cutoff;
    SpeechEEGCoherence_Results.lowpass_cutoff = lowpass_cutoff;
    SpeechEEGCoherence_Results.RandomSeed = options.RandomSeed;
    SpeechEEGCoherence_Results.label = EEG_Speech_prepped.label;
    SpeechEEGCoherence_Results.Chanlocs = EEG_Speech_prepped.Chanlocs;
    SpeechEEGCoherence_Results.Condition = EEG_Speech_prepped.Condition;
    SpeechEEGCoherence_Results.Subj_id = EEG_Speech_prepped.Subj_id;
    SpeechEEGCoherence_Results.audio_file = EEG_Speech_prepped.audio_file;
 
    % if we NumShuffles is 0 or negative, skip the shuffle test and return only the true band-averaged MSC.
    if options.NumShuffles < 1
        fprintf(' -- NumShuffles < 1, skipping shuffle test and returning only true band-averaged MSC.\n');
        return;     % break out of the function early, returning only the true band-averaged MSC and metadata.
    end

    % ------------------------------------------------------------------------------------------------------------
    % Step 3: Shuffle loop. Each iteration draws a fresh random derangement of the trial order,
    % computes trial-averaged CSD/PSD for that mismatched pairing,
    % and band-averages the result via BandAverageMSC.m.
    % ------------------------------------------------------------------------------------------------------------
    exceed_count = zeros(Num_channels, 1);
    null_sum = zeros(Num_channels, 1);
    null_sumsq = zeros(Num_channels, 1);
    PairingsUsed = zeros(options.NumShuffles, Num_trials);
    Null_MSC = zeros(options.NumShuffles, Num_channels);

    progress_step = max(1, round(options.NumShuffles / 10));
    for s = 1:options.NumShuffles
        pairing_s = derange(1:Num_trials);
        PairingsUsed(s, :) = pairing_s;

        NullResults_s = ComputeCoherenceForPairing(EEG_Speech_prepped, pairing_s);
        NullBandResults_s = BandAverageMSC(NullResults_s, band_low, band_high, AvgMethod=options.AvgMethod);
        msc_band_s = NullBandResults_s.msc_band; % [Num_channels x 1]

        exceed_count = exceed_count + double(msc_band_s >= True_MSC);
        null_sum = null_sum + msc_band_s;
        null_sumsq = null_sumsq + msc_band_s.^2;
        Null_MSC(s, :) = msc_band_s(:)';

        if mod(s, progress_step) == 0
            fprintf('   shuffle %d / %d complete\n', s, options.NumShuffles);
        end
    end

    % ------------------------------------------------------------------------------------------------------------
    % Step 4: summarize the null distribution and compute p-values / z-scores.
    % ------------------------------------------------------------------------------------------------------------
    NumShuffles = options.NumShuffles;
    Null_Mean = null_sum / NumShuffles;
    Null_Var = (null_sumsq - NumShuffles * Null_Mean.^2) / max(NumShuffles - 1, 1);
    Null_Std = sqrt(max(Null_Var, 0)); % guard against tiny negative values from floating-point roundoff

    PValues = (1 + exceed_count) / (NumShuffles + 1);
    ZScores = (True_MSC - Null_Mean) ./ Null_Std;



    % ------------------------------------------------------------------------------------------------------------
    % Step 5: package results.
    % ------------------------------------------------------------------------------------------------------------
    SpeechEEGCoherence_Results.PValues = PValues; % [Num_channels x 1] one-sided permutation p-values, per channel.
    SpeechEEGCoherence_Results.ZScores = ZScores;   % [Num_channels x 1] (True_MSC - Null_Mean) ./ Null_Std, per channel (a parametric-flavored companion to the p-values -- useful for ranking/plotting, but the p-values are the actual significance test).
    SpeechEEGCoherence_Results.Null_Mean = Null_Mean;   % [Num_channels x 1] mean of the null band-averaged MSC distribution.
    SpeechEEGCoherence_Results.Null_Std = Null_Std;  % [Num_channels x 1] std of the null band-averaged MSC distribution.
    SpeechEEGCoherence_Results.Null_MSC = Null_MSC;  % [NumShuffles x Num_channels] the full null distribution (band-averaged MSC is a single number per channel per shuffle, so this is small -- a few hundred KB even for thousands of shuffles -- and is always returned).
    SpeechEEGCoherence_Results.PairingsUsed = PairingsUsed;   % [NumShuffles x Num_trials] the derangement used for each shuffle.
    SpeechEEGCoherence_Results.NumShuffles = NumShuffles;  % number of random derangements drawn for the null distribution.
    SpeechEEGCoherence_Results.TestDescription = sprintf(...
        ['One-sided permutation test on band-averaged [%.2f %.2f] Hz MSC (AvgMethod=%s): ' ...
         'p = P(null band-MSC >= true band-MSC), estimated from %d random derangements of the ' ...
         'trial order (no trial paired with its own speech epoch). ' ...
         'p = (1 + #shuffles with null >= true) / (NumShuffles + 1). Not corrected for multiple comparisons across channels.'], ...
        band_low, band_high, options.AvgMethod, NumShuffles);
    fprintf('=== ShuffleTest_SpeechEEGCoherence complete for subject %s ===\n', eeg_struct.Subj_id);
end


SpeechEEGCoherence_Results.PValues = zeros(Num_channels, 1); % [Num_channels x 1] one-sided permutation p-values, per channel.
SpeechEEGCoherence_Results.ZScores = zeros(Num_channels, 1);   % [Num_channels x 1] (True_MSC - Null_Mean) ./ Null_Std, per channel (a parametric-flavored companion to the p-values -- useful for ranking/plotting, but the p-values are the actual significance test).
SpeechEEGCoherence_Results.Null_Mean = zeros(Num_channels, 1);   % [Num_channels x 1] mean of the null band-averaged MSC distribution.
SpeechEEGCoherence_Results.Null_Std = zeros(Num_channels, 1);  % [Num_channels x 1] std of the null band-averaged MSC distribution.
% SpeechEEGCoherence_Results.Null_MSC = zeros(NumShuffles, Num_channels);  % [NumShuffles x Num_channels] the full null distribution (band-averaged MSC is a single number per channel per shuffle, so this is small -- a few hundred KB even for thousands of shuffles -- and is always returned).
% SpeechEEGCoherence_Results.PairingsUsed = PairingsUsed;   % [NumShuffles x Num_trials] the derangement used for each shuffle.
% SpeechEEGCoherence_Results.NumShuffles = NumShuffles;  % number of random derangements drawn for the null distribution.
SpeechEEGCoherence_Results.TestDescription = sprintf(...
