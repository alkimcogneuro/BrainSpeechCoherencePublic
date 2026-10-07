function [EEG_Speech_prepped] = PreprocessSpeechEEGForCoherence(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % PREPROCESS SPEECH + EEG DATA FOR COHERENCE ANALYSIS (pairing-independent step)
    % ===============================================================================================================
    % This function preprocesses the EEG and speech data for coherence analysis, 
    % including z-scoring, band-pass filtering, and extracting the speech envelope. 
    % 
    % It is designed to be run once per subject. 
    % The results can be reused for multiple trial pairings, as in a shuffle test, 
    % which repeatedly pairs each EEG trial with the speech from a different trial,
    % as in ShuffleTest_SpeechEEGCoherence.m. 
    % This allows us to avoid repeating the expensive preprocessing steps.
    % The code will do the following:
    %   - z-score the EEG data (per channel, across all samples/trials)
    %   - z-score the speech amplitude waveform
    %   - band-pass filtering each trial's EEG epoch, for each channel (filtfilt -- the slow step)
    %   - extracting each trial's OWN speech epoch (at that trial's OWN original onset latency),
    %     then computing its Hilbert envelope, band-pass filtering it, and downsampling it to the
    %     EEG sampling rate (hilbert + filtfilt + resample -- also slow)
    %

    % "Pairing-independent" means: nothing here depends on which trial's speech epoch ends up
    % analyzed against which trial's EEG epoch. That includes:
    % 
    % None of that changes if we later decide to analyze EEG trial i against a *different* trial's
    % speech epoch (which is exactly what a shuffle/permutation test does). So we do it once here,
    % and hand the results to ComputeCoherenceForPairing.m, which is cheap to call repeatedly.
    %
    % IMPORTANT -- keep this in sync with AnalyzeSpeechEEGCoherence_Dataset.m:
    % This function intentionally duplicates the validation / epoching / normalization logic at
    % the top of AnalyzeSpeechEEGCoherence_Dataset.m, rather than modifying that (already-validated,
    % already-in-use) function. If you change how epoching, z-scoring, or filtering works in one of
    % these two files, make the same change in the other, or the shuffle test will silently stop
    % matching the main analysis pipeline.
    %
    % ARGUMENTS: 
    %   same meaning as the corresponding arguments in AnalyzeSpeechEEGCoherence_Dataset.m.
    %   eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options.StartTimeOffset, options.EpochDuration
    %
    % Welch-windowing arguments 
    % CSD/PSD are estimated per trial via WelchCrossSpectralDensity.m:  
    % several overlapping windows within the trial, averaged
    % together, rather than one window spanning the whole trial; 
    % see that file's header, and
    % RecommendWelchWindowLength.m, for the full rationale):
    %   options.WindowDurationSec (default NaN): analysis window length, in seconds, used for
    %       every trial's CSD/PSD estimate. NaN (the default) means "use the floor window length
    %       needed to resolve options.MinFreqHz at options.MinCyclesPerWindow cycles" -- i.e., the
    %       shortest window that still resolves your lowest frequency of interest, which maximizes
    %       the number of segments/trial. Pass this explicitly once you've decided on a value for a
    %       given dataset (e.g., via RecommendWelchWindowLength.m) rather than relying on the auto
    %       floor for every dataset -- trial length varies a lot across datasets (this codebase has
    %       used both ~30 sec and ~4 sec trials), and the right window length for a 30 sec trial
    %       (plenty of room for a longer window with many segments) is not the right window length
    %       for a 4 sec trial (already close to the frequency-resolution floor).
    %   options.OverlapFraction (default 0.5): fraction of window overlap between consecutive
    %       segments within a trial. Must be in [0, 0.5] -- see WelchCrossSpectralDensity.m header
    %       for why more than 50%% overlap isn't accepted (segments become too correlated with their
    %       neighbor to behave like separate observations in the average).
    %   options.MinFreqHz (default = highpass_cutoff): the lowest frequency this analysis needs to
    %       resolve reliably, used only to derive the floor window length when WindowDurationSec is
    %       NaN. Defaults to highpass_cutoff since there's no point resolving frequencies you've
    %       already filtered out.
    %   options.MinCyclesPerWindow (default 3): minimum full cycles of MinFreqHz required inside one
    %       analysis window; see RecommendWelchWindowLength.m.
    %
    % Note: there is no RandomizeOnsets option here. In the original pipeline, RandomizeOnsets
    % worked by deranging *which onset latency* is used to extract each trial's speech epoch --
    % but since every trial's speech epoch has the same fixed duration for a given subject,
    % extracting trial i's speech using trial pi(i)'s onset latency produces EXACTLY the same
    % samples as trial pi(i)'s own (originally-extracted) speech epoch. So instead of re-extracting
    % per shuffle, we extract each trial's own speech epoch once here (columns of
    % speech_envelope_trials, indexed by each trial's OWN original onset), and shuffling is done
    % downstream by permuting which column of speech_envelope_trials gets paired with which EEG
    % trial (see ComputeCoherenceForPairing.m / ShuffleTest_SpeechEEGCoherence.m). This is
    % mathematically identical to the original RandomizeOnsets approach, just reorganized so the
    % expensive per-epoch work (filtering, Hilbert transform, resampling) happens once instead of
    % once per shuffle.
    %
    % Returns a structure EEG_Speech_prepped with fields:
    %   .eeg_epoch_bpf          [Num_channels x Num_samples_eeg x Num_trials] band-pass filtered, z-scored EEG
    %   .speech_envelope_trials [Num_samples_eeg x Num_trials] preprocessed speech envelope for each
    %                           trial's OWN original speech epoch (column j = trial j's own epoch)
    %   .nfft, .num_freqs, .Freqs   FFT / frequency-axis parameters -- as of 2026-09-17, .nfft is
    %                           the WELCH ANALYSIS WINDOW length in samples (see WindowDurationSec
    %                           above), NOT the full trial length -- .num_freqs/.Freqs are derived
    %                           from that window length, since that's what actually sets the
    %                           frequency resolution now that CSD/PSD come from averaging multiple
    %                           windows per trial rather than one window spanning the whole trial.
    %   .WelchWindowSamples     same value as .nfft -- explicit alias, kept for readability at call
    %                           sites that specifically care about "the Welch window length" rather
    %                           than "the nfft used," even though today they're the same number.
    %   .WelchOverlapFraction   the overlap fraction actually used (options.OverlapFraction, echoed
    %                           back for bookkeeping/provenance).
    %   .WelchNumSegmentsPerTrial   the number of Welch segments each trial's CSD/PSD estimate
    %                           actually averaged over, given WelchWindowSamples/WelchOverlapFraction
    %                           and this dataset's trial length (see RecommendWelchWindowLength.m).
    %   .eeg_Fs, .speech_Fs_downsampled
    %   .Num_channels, .Num_trials
    %   .highpass_cutoff, .lowpass_cutoff  (passthrough, for bookkeeping)
    %   .label, .Chanlocs, .Condition, .Subj_id, .audio_file   (passthrough metadata, same source
    %                           fields used by AnalyzeSpeechEEGCoherence_Dataset's CoherenceResults)
    %   .OnsetLatency_used      the (possibly StartTimeOffset-shifted) onset latencies actually used
    %                           to extract each trial's own speech epoch -- kept for bookkeeping/debugging.
    % ===============================================================================================================
    arguments
        eeg_struct
        speech_rawdata
        highpass_cutoff (1,1) double {mustBeReal} = 2
        lowpass_cutoff  (1,1) double {mustBeReal} = 35
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN
        options.EpochDuration (1,1) double {mustBeReal} = NaN
        options.WindowDurationSec (1,1) double {mustBeReal} = NaN
        options.OverlapFraction (1,1) double {mustBeReal} = 0.5
        options.MinFreqHz (1,1) double {mustBeReal} = highpass_cutoff
        options.MinCyclesPerWindow (1,1) double {mustBeReal} = 3
    end

    fprintf(' -- PreprocessSpeechEEGForCoherence() for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Speech file is %s\n', eeg_struct.audio_file);
    fprintf(' -- Filter parameters: HighPass = %.2f Hz, LowPass = %.2f Hz\n', highpass_cutoff, lowpass_cutoff);
    fprintf(' -- EEG sampling rate: %d Hz, Speech sampling rate: %d Hz\n', eeg_struct.Fs, speech_rawdata.Fs);
    fprintf(' -- Number of EEG trials: %d, Number of EEG channels: %d\n', eeg_struct.Num_trials, eeg_struct.Num_channels);

    % --------------------------------------------------------------------------------------------------------------
    % Check the epoch timing parameters. 
    % The rationale for these checks is provided in more detail in AnalyzeSpeechEEGCoherence_Dataset.m 
    % --------------------------------------------------------------------------------------------------------------
    if ~isnan(options.StartTimeOffset) || ~isnan(options.EpochDuration)
        if ~isnan(options.StartTimeOffset) && isnan(options.EpochDuration)
            error('If StartTimeOffset is provided, EpochDuration must also be provided.');
        end
        if ~isnan(options.EpochDuration) && isnan(options.StartTimeOffset)
            error('If EpochDuration is provided, StartTimeOffset must also be provided.');
        end
        if ~isnan(options.StartTimeOffset) && (options.StartTimeOffset < 0)
            error('The specified StartTimeOffset of %.2f seconds is negative.', options.StartTimeOffset);
        end
        if ~isnan(options.EpochDuration) && (options.EpochDuration > (size(eeg_struct.Data, 2) / eeg_struct.Fs))
            error('The specified EpochDuration of %.2f seconds is too long for the length of the EEG epochs (%.2f seconds).', ...
                options.EpochDuration, size(eeg_struct.Data, 2) / eeg_struct.Fs);
        end
        if ~isnan(options.StartTimeOffset) && ~isnan(options.EpochDuration) && (options.StartTimeOffset > options.EpochDuration)
            error('StartTimeOffset (%.2f s) is greater than EpochDuration (%.2f s).', options.StartTimeOffset, options.EpochDuration);
        end
        start_offset_samples = round(options.StartTimeOffset * eeg_struct.Fs);
        epoch_duration_samples = round(options.EpochDuration * eeg_struct.Fs);
        sample_range_restricted = start_offset_samples:(start_offset_samples + epoch_duration_samples - 1);
        eeg_struct.Data = eeg_struct.Data(:, sample_range_restricted, :);
        eeg_struct.OnsetLatency = eeg_struct.OnsetLatency + options.StartTimeOffset;
    end

    if any(isnan(eeg_struct.OnsetLatency))
        error('EEG_struct.OnsetLatency contains NaN values. Please check the EEG data structure.');
    end

    Num_channels = size(eeg_struct.Data, 1);
    Num_samples_eeg = size(eeg_struct.Data, 2);
    Num_trials = size(eeg_struct.Data, 3);

    if Num_channels > 64
        fprintf(' -- EEG data has %d channels; restricting to the first 64 to match the standard montage.\n', Num_channels);
        eeg_struct.Data = eeg_struct.Data(1:64, :, :);
        eeg_struct.Chanlocs = eeg_struct.Chanlocs(1:64);
        Num_channels = 64;
    end

    % ------------------------------------------------------------------------------------------------------------
    % Z-score EEG (per channel, across all samples and trials) and speech (across the whole recording).
    % Identical logic to AnalyzeSpeechEEGCoherence_Dataset.m.
    % ------------------------------------------------------------------------------------------------------------
    eegdata_2d = reshape(eeg_struct.Data, Num_channels, []);  % this call to reshape will not copy the data. 
    % It just creates a new view of the same underlying array. So this is memory-efficient.
    % The resulting eegdata_2d is [Num_channels x (Num_samples_eeg * Num_trials)]
    % that will allow us to compute the mean and std for each channel across all samples and trials, 
    % and then z-score each channel's data.
    % Then we'll reshape it back to the original [Num_channels x Num_samples_eeg x Num_trials] format
    % for further processing.

    mu    = mean(eegdata_2d, 2);    % mean across the second dimension (samples and trials) for each channel
    sigma = std(eegdata_2d, 0, 2);  % standard deviation across the second dimension (samples and trials) for each channel
    eeg_data_zscored = (eegdata_2d - mu) ./ sigma; % z-score each channel's data
    eeg_struct.Data = reshape(eeg_data_zscored, Num_channels, Num_samples_eeg, Num_trials);        % reshape back to the original format
    % and replace the original data with the z-scored version


    speech_amplitudes = speech_rawdata.Amplitudes;
    speech_amplitudes_zscored = (speech_amplitudes - mean(speech_amplitudes)) / std(speech_amplitudes); % z-score the speech amplitudes
    speech_rawdata.Amplitudes = speech_amplitudes_zscored;  % replace the original amplitudes with the z-scored version
    num_speech_samples = length(speech_rawdata.Amplitudes);

    % --------------------------------------------------------------------------------------------------------------
    % Welch-style windowing parameters: determine (and validate) the analysis window length used
    % for CSD/PSD estimation. Each trial's CSD/PSD is now estimated from multiple overlapping
    % windows within the trial (see WelchCrossSpectralDensity.m) rather than one window spanning
    % the whole trial. Validated ONCE here (not per-trial/per-shuffle) so a bad configuration fails
    % fast with a clear message instead of failing deep inside a long shuffle loop.
    % --------------------------------------------------------------------------------------------------------------
    trial_duration_sec = Num_samples_eeg / eeg_struct.Fs;
    if isnan(options.WindowDurationSec)
        WelchWindowInfo = RecommendWelchWindowLength(trial_duration_sec, eeg_struct.Fs, options.MinFreqHz, ...
            MinCyclesPerWindow=options.MinCyclesPerWindow, OverlapFraction=options.OverlapFraction);
        fprintf([' -- options.WindowDurationSec not specified; using the floor window length for ' ...
                 '%.3g Hz (%g cycles/window): %.3f sec (%d samples), giving %d segment(s)/trial.\n'], ...
                 options.MinFreqHz, options.MinCyclesPerWindow, WelchWindowInfo.WindowDurationSec, ...
                 WelchWindowInfo.WindowSamples, WelchWindowInfo.NumSegments);
    else
        WelchWindowInfo = RecommendWelchWindowLength(trial_duration_sec, eeg_struct.Fs, options.MinFreqHz, ...
            MinCyclesPerWindow=options.MinCyclesPerWindow, OverlapFraction=options.OverlapFraction, ...
            WindowDurationSec=options.WindowDurationSec);
        fprintf(' -- Welch analysis window: %.3f sec (%d samples), %.0f%% overlap, %d segment(s)/trial.\n', ...
                 WelchWindowInfo.WindowDurationSec, WelchWindowInfo.WindowSamples, options.OverlapFraction*100, ...
                 WelchWindowInfo.NumSegments);
    end

    nfft = WelchWindowInfo.WindowSamples;   % the WELCH WINDOW/segment length now, not the full trial length
    if mod(nfft, 2) == 0
        num_freqs = nfft/2 + 1;
    else
        num_freqs = (nfft + 1)/2;
    end
    all_freqs = (0:nfft-1) * (eeg_struct.Fs / nfft);
    Freqs = all_freqs(1:num_freqs);

    speech_epoch_duration_seconds = Num_samples_eeg / eeg_struct.Fs;
    speech_epoch_duration_samples = round(speech_epoch_duration_seconds * speech_rawdata.Fs);

    % ------------------------------------------------------------------------------------------------------------
    % Preallocate outputs.
    % ------------------------------------------------------------------------------------------------------------
    eeg_epoch_bpf = zeros(Num_channels, Num_samples_eeg, Num_trials);
    speech_envelope_trials = zeros(Num_samples_eeg, Num_trials);
    OnsetLatency_used = eeg_struct.OnsetLatency(:)';    % make sure it's a row vector for consistent indexing later

    % ------------------------------------------------------------------------------------------------------------
    % Loop over trials: extract & preprocess each trial's OWN speech epoch, and band-pass filter each
    % trial's EEG epoch for every channel. This is the expensive step (filtfilt, hilbert, resample) --
    % we do it exactly once per trial, regardless of how many shuffles will later reuse the result.
    % ------------------------------------------------------------------------------------------------------------
    for trial_idx = 1:Num_trials
        speech_onset_latency = OnsetLatency_used(trial_idx);
        if speech_onset_latency == 0
            speech_onset_idx = 1;
        else
            % translate the speech onset latency (in seconds) to an index in the speech amplitude array (in samples)
            speech_onset_idx = round(speech_onset_latency * speech_rawdata.Fs);
        end
        speech_offset_idx = speech_onset_idx + speech_epoch_duration_samples - 1;

        if speech_onset_idx < 1
            error('PreprocessSpeechEEGForCoherence:OnsetBeforeRecordingStart', ...
                'Trial %d: computed speech onset index (%d samples, %.3f sec) is before the start of the speech recording.', ...
                trial_idx, speech_onset_idx, speech_onset_latency);
        end
        if speech_offset_idx > num_speech_samples
            error('PreprocessSpeechEEGForCoherence:OnsetTooCloseToRecordingEnd', ...
                ['Trial %d: speech epoch requires samples %d:%d, but the speech recording is only %d samples long. ' ...
                 'Onset latency = %.3f sec, required epoch duration = %.3f sec (%d samples).'], ...
                trial_idx, speech_onset_idx, speech_offset_idx, num_speech_samples, ...
                speech_onset_latency, speech_epoch_duration_samples / speech_rawdata.Fs, speech_epoch_duration_samples);
        end

        speech_epoch = speech_rawdata.Amplitudes(speech_onset_idx:speech_offset_idx);

        % returns the amplitude envelope, downsampled to match the EEG sampling rate, 
        % as well as the instantaneous phase values of the speech envelope (and reversed speech envelope).
        % currenlty, we only use the envelope, but the phase values are returned for potential future use.
        [Speech_Struct] = preprocess_speech_epoch(speech_epoch, speech_rawdata.Fs, ...
                                                  highpass_cutoff, lowpass_cutoff, Num_samples_eeg);
        speech_envelope_trials(:, trial_idx) = Speech_Struct.envelope(:);

        % iterate through each channel and band-pass filter the EEG epoch for that channel and trial
        for ch_idx = 1:Num_channels
            eeg_epoch = squeeze(eeg_struct.Data(ch_idx, :, trial_idx)); % one channel's data for one trial
            eeg_epoch_bpf(ch_idx, :, trial_idx) = band_pass_filt(eeg_epoch, eeg_struct.Fs, highpass_cutoff, lowpass_cutoff);
        end
        fprintf('   preprocessed trial %d / %d\n', trial_idx, Num_trials);
    end

    % ------------------------------------------------------------------------------------------------------------
    % Package results.
    % ------------------------------------------------------------------------------------------------------------
    EEG_Speech_prepped = struct();
    EEG_Speech_prepped.eeg_epoch_bpf = eeg_epoch_bpf;
    EEG_Speech_prepped.speech_envelope_trials = speech_envelope_trials;
    EEG_Speech_prepped.nfft = nfft;
    EEG_Speech_prepped.num_freqs = num_freqs;
    EEG_Speech_prepped.Freqs = Freqs;
    EEG_Speech_prepped.WelchWindowSamples = nfft;               % same value as .nfft -- explicit alias for readability at call sites
    EEG_Speech_prepped.WelchOverlapFraction = options.OverlapFraction;
    EEG_Speech_prepped.WelchNumSegmentsPerTrial = WelchWindowInfo.NumSegments;
    EEG_Speech_prepped.eeg_Fs = eeg_struct.Fs;
    EEG_Speech_prepped.speech_Fs_downsampled = Speech_Struct.fs; % all trials share the same downsampled rate (== eeg_Fs, by construction)
    EEG_Speech_prepped.Num_channels = Num_channels;
    EEG_Speech_prepped.Num_trials = Num_trials;
    EEG_Speech_prepped.highpass_cutoff = highpass_cutoff;
    EEG_Speech_prepped.lowpass_cutoff = lowpass_cutoff;
    EEG_Speech_prepped.label = eeg_struct.Subj_id;
    EEG_Speech_prepped.Chanlocs = eeg_struct.Chanlocs;
    EEG_Speech_prepped.Condition = eeg_struct.Condition;
    EEG_Speech_prepped.Subj_id = eeg_struct.Subj_id;
    EEG_Speech_prepped.audio_file = eeg_struct.audio_file;
    EEG_Speech_prepped.OnsetLatency_used = OnsetLatency_used;
end
