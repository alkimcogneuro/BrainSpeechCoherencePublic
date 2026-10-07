function [ShuffleTestResults] = ShuffleTest_SpeechEEGCoherence(eeg_struct, speech_rawdata, band_low, band_high, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % SHUFFLE-BASED SIGNIFICANCE TEST FOR BAND-AVERAGED SPEECH-EEG MAGNITUDE SQUARED COHERENCE
    % ===============================================================================================================
    % Tests, for one subject, at every EEG channel, whether the TRUE (correctly-paired) trial-averaged,
    % BAND-AVERAGED magnitude squared coherence (MSC) between the speech envelope and the EEG is
    % greater than expected by chance -- i.e. exactly one permutation test per channel, run on the
    % same band-averaged quantity that BandAverageMSC.m produces, rather than one test per channel
    % per frequency bin.
    %    %
    % "Expected by chance" is defined as follows: 
    % Take a random derangement of the trial order (derange.m -- a permutation with
    % no trial mapped to itself), so every EEG trial is analyzed against a DIFFERENT trial's speech
    % epoch, band-average the resulting MSC (via BandAverageMSC.m), and repeat. Doing this many times
    % builds a null distribution of band-averaged MSC values, per channel, against which the true
    % band-averaged MSC is compared.
    %
    % This function calls PreprocessSpeechEEGForCoherence.m ONCE (the expensive step: z-scoring,
    % band-pass filtering every trial/channel, and extracting + filtering every trial's speech
    % envelope), then calls ComputeCoherenceForPairing.m (cheap: FFT + averaging, still computed at
    % full frequency resolution -- BandAverageMSC needs the full CSD/PSD spectra to average
    % correctly, see its header) once for the true pairing and once per shuffle, band-averaging each
    % result. See PreprocessSpeechEEGForCoherence.m / ComputeCoherenceForPairing.m for the full
    % rationale on why preprocessing is separated from the pairing-dependent computation.
    %
    % Significance is one-sided: p(channel) = P(null band-MSC >= true band-MSC), using the standard
    % permutation-test estimator that never returns exactly zero:
    %     p = (1 + #{shuffles with null band-MSC >= true band-MSC}) / (NumShuffles + 1)
    % (Phipson & Smyth, 2010.) The smallest possible p-value is therefore 1/(NumShuffles+1), so
    % NumShuffles should be chosen with that resolution in mind (e.g. 1000 shuffles -> smallest
    % possible p is ~0.001).
    %
    % NOTE on multiple comparisons: this test is still run independently at every channel (e.g. 64
    % tests for a 64-channel montage), so we might consider FDR correction (e.g. Benjamini-Hochberg) across
    % channels before declaring significance, even though the frequency dimension has been collapsed.
    %
    % Arguments:
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
    %   options.NumShuffles (default 1000) : number of random derangements to draw for the null
    %     distribution.
    %   options.RandomSeed (default NaN = don't seed) : if provided, calls rng(RandomSeed) before
    %     drawing shuffles, for reproducibility. NOTE: this reseeds MATLAB's global random stream.
    %   options.AvgMethod (default 'power-weighted') : passed through to BandAverageMSC.m -- either
    %     'power-weighted' (recommended, matches how CSD/PSD are already averaged across trials in
    %     this pipeline) or 'simple'. See BandAverageMSC.m for the rationale.
    %
    % Returns ShuffleTestResults, a struct with fields:
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
        options.NumShuffles (1,1) double {mustBePositive, mustBeInteger} = 1000
        options.RandomSeed (1,1) double {mustBeReal} = NaN
        options.AvgMethod (1,:) char = 'power-weighted'
    end

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
    % Step 1: expensive, pairing-independent preprocessing -- done once. Still computed at full
    % frequency resolution inside ComputeCoherenceForPairing below; only the accumulation/statistics
    % in this file are now band-level.
    % ------------------------------------------------------------------------------------------------------------
    PreppedData = PreprocessSpeechEEGForCoherence(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, ...
        StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration);

    Num_channels = PreppedData.Num_channels;
    Num_trials = PreppedData.Num_trials;

    % ------------------------------------------------------------------------------------------------------------
    % Step 2: the TRUE / observed pairing (trial i's own EEG paired with trial i's own speech),
    % band-averaged.
    % ------------------------------------------------------------------------------------------------------------
    fprintf(' -- Computing true (observed) band-averaged coherence...\n');
    TrueResults = ComputeCoherenceForPairing(PreppedData, 1:Num_trials);
    [True_MSC, band_freqs] = BandAverageMSC(TrueResults, band_low, band_high, AvgMethod=options.AvgMethod);

    % ------------------------------------------------------------------------------------------------------------
    % Step 3: shuffle loop. Each iteration draws a fresh random derangement of the trial order
    % (matching the original RandomizeOnsets logic exactly -- see PreprocessSpeechEEGForCoherence.m
    % header for the equivalence argument), recomputes trial-averaged CSD/PSD for that mismatched
    % pairing, and band-averages the result via BandAverageMSC.m.
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

        NullResults_s = ComputeCoherenceForPairing(PreppedData, pairing_s);
        msc_band_s = BandAverageMSC(NullResults_s, band_low, band_high, AvgMethod=options.AvgMethod); % [Num_channels x 1]

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
    N = options.NumShuffles;
    Null_Mean = null_sum / N;
    Null_Var = (null_sumsq - N * Null_Mean.^2) / max(N - 1, 1);
    Null_Std = sqrt(max(Null_Var, 0)); % guard against tiny negative values from floating-point roundoff

    PValues = (1 + exceed_count) / (N + 1);
    ZScores = (True_MSC - Null_Mean) ./ Null_Std;

    % ------------------------------------------------------------------------------------------------------------
    % Step 5: package results.
    % ------------------------------------------------------------------------------------------------------------
    ShuffleTestResults = struct();
    ShuffleTestResults.True_MSC = True_MSC;
    ShuffleTestResults.PValues = PValues;
    ShuffleTestResults.ZScores = ZScores;
    ShuffleTestResults.Null_Mean = Null_Mean;
    ShuffleTestResults.Null_Std = Null_Std;
    ShuffleTestResults.Null_MSC = Null_MSC;
    ShuffleTestResults.PairingsUsed = PairingsUsed;
    ShuffleTestResults.band_low = band_low;
    ShuffleTestResults.band_high = band_high;
    ShuffleTestResults.band_freqs = band_freqs;
    ShuffleTestResults.AvgMethod = options.AvgMethod;
    ShuffleTestResults.NumShuffles = N;
    ShuffleTestResults.nfft = PreppedData.nfft;
    ShuffleTestResults.eeg_Fs = PreppedData.eeg_Fs;
    ShuffleTestResults.highpass_cutoff = highpass_cutoff;
    ShuffleTestResults.lowpass_cutoff = lowpass_cutoff;
    ShuffleTestResults.RandomSeed = options.RandomSeed;
    ShuffleTestResults.label = PreppedData.label;
    ShuffleTestResults.Chanlocs = PreppedData.Chanlocs;
    ShuffleTestResults.Condition = PreppedData.Condition;
    ShuffleTestResults.Subj_id = PreppedData.Subj_id;
    ShuffleTestResults.audio_file = PreppedData.audio_file;
    ShuffleTestResults.TestDescription = sprintf(...
        ['One-sided permutation test on band-averaged [%.2f %.2f] Hz MSC (AvgMethod=%s): ' ...
         'p = P(null band-MSC >= true band-MSC), estimated from %d random derangements of the ' ...
         'trial order (no trial paired with its own speech epoch). ' ...
         'p = (1 + #shuffles with null >= true) / (NumShuffles + 1). Not corrected for multiple comparisons across channels.'], ...
        band_low, band_high, options.AvgMethod, N);

    fprintf('=== ShuffleTest_SpeechEEGCoherence complete for subject %s ===\n', eeg_struct.Subj_id);
end
