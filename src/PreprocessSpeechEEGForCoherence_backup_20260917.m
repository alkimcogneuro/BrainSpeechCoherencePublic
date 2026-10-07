function [PreppedData] = PreprocessSpeechEEGForCoherence(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % PREPROCESS SPEECH + EEG DATA FOR COHERENCE ANALYSIS (pairing-independent step)
    % ===============================================================================================================
    % This function factors the *expensive, pairing-independent* preprocessing steps out of
    % AnalyzeSpeechEEGCoherence_Dataset.m, so they can be run ONCE per subject and then reused
    % many times by ShuffleTest_SpeechEEGCoherence.m (once for the true/observed trial pairing,
    % and once per shuffle).
    %
    % "Pairing-independent" means: nothing here depends on which trial's speech epoch ends up
    % analyzed against which trial's EEG epoch. That includes:
    %   - z-scoring the EEG data (per channel, across all samples/trials)
    %   - z-scoring the speech amplitude waveform
    %   - band-pass filtering each trial's EEG epoch, for each channel (filtfilt -- the slow step)
    %   - extracting each trial's OWN speech epoch (at that trial's OWN original onset latency),
    %     then computing its Hilbert envelope, band-pass filtering it, and downsampling it to the
    %     EEG sampling rate (hilbert + filtfilt + resample -- also slow)
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
    % Arguments: same meaning as the corresponding arguments in AnalyzeSpeechEEGCoherence_Dataset.m.
    %   eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options.StartTimeOffset, options.EpochDuration
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
    % Returns a structure PreppedData with fields:
    %   .eeg_epoch_bpf          [Num_channels x Num_samples_eeg x Num_trials] band-pass filtered, z-scored EEG
    %   .speech_envelope_trials [Num_samples_eeg x Num_trials] preprocessed speech envelope for each
    %                           trial's OWN original speech epoch (column j = trial j's own epoch)
    %   .nfft, .num_freqs, .Freqs   FFT / frequency-axis parameters (identical meaning to
    %                           AnalyzeSpeechEEGCoherence_Dataset's CoherenceResults.nfft / .Freqs)
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
    end

    fprintf(' -- PreprocessSpeechEEGForCoherence() for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Speech file is %s\n', eeg_struct.audio_file);
    fprintf(' -- Filter parameters: HighPass = %.2f Hz, LowPass = %.2f Hz\n', highpass_cutoff, lowpass_cutoff);
    fprintf(' -- EEG sampling rate: %d Hz, Speech sampling rate: %d Hz\n', eeg_struct.Fs, speech_rawdata.Fs);
    fprintf(' -- Number of EEG trials: %d, Number of EEG channels: %d\n', eeg_struct.Num_trials, eeg_struct.Num_channels);

    % --------------------------------------------------------------------------------------------------------------
    % Epoch-window validation and application (StartTimeOffset / EpochDuration).
    % Copied verbatim (in spirit) from AnalyzeSpeechEEGCoherence_Dataset.m -- see that file for the
    % detailed rationale on each check.
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
    eegdata_2d = reshape(eeg_struct.Data, Num_channels, []);  % this call to reshape will not copy the data, 
    % it just creates a new view of the same underlying array. So this is memory-efficient.
    % the resulting eegdata_2d is [Num_channels x (Num_samples_eeg * Num_trials)]
    % that will allow us to compute the mean and std for each channel across all samples and trials, 
    % and then z-score each channel's data.
    % then we'll reshape it back to [Num_channels x Num_samples_eeg x Num_trials] for further processing.

    mu    = mean(eegdata_2d, 2);
    sigma = std(eegdata_2d, 0, 2);
    eeg_data_norm = (eegdata_2d - mu) ./ sigma;
    eeg_struct.Data = reshape(eeg_data_norm, Num_channels, Num_samples_eeg, Num_trials);

    speech_amplitudes = speech_rawdata.Amplitudes;
    speech_amplitudes_norm = (speech_amplitudes - mean(speech_amplitudes)) / std(speech_amplitudes);
    speech_rawdata.Amplitudes = speech_amplitudes_norm;
    num_speech_samples = length(speech_rawdata.Amplitudes);

    nfft = Num_samples_eeg;
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
    OnsetLatency_used = eeg_struct.OnsetLatency(:)';

    % ------------------------------------------------------------------------------------------------------------
    % Loop over trials: extract + preprocess each trial's OWN speech epoch, and band-pass filter each
    % trial's EEG epoch for every channel. This is the expensive step (filtfilt, hilbert, resample) --
    % we do it exactly once per trial, regardless of how many shuffles will later reuse the result.
    % ------------------------------------------------------------------------------------------------------------
    for trial_idx = 1:Num_trials
        speech_onset_latency = OnsetLatency_used(trial_idx);
        if speech_onset_latency == 0
            speech_onset_idx = 1;
        else
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
        % and the instantaneous phase values of the speech envelope (and reversed speech envelope).

        [Speech_Struct] = preprocess_speech_epoch(speech_epoch, speech_rawdata.Fs, ...
                                                  highpass_cutoff, lowpass_cutoff, Num_samples_eeg);
        speech_envelope_trials(:, trial_idx) = Speech_Struct.envelope(:);

        for ch_idx = 1:Num_channels
            eeg_epoch = squeeze(eeg_struct.Data(ch_idx, :, trial_idx)); % one channel's data for one trial
            if size(eeg_epoch, 2) ~= 1
                eeg_epoch = eeg_epoch';
            end
            eeg_epoch_bpf(ch_idx, :, trial_idx) = band_pass_filt(eeg_epoch, eeg_struct.Fs, highpass_cutoff, lowpass_cutoff);
        end
        fprintf('   preprocessed trial %d / %d\n', trial_idx, Num_trials);
    end

    % ------------------------------------------------------------------------------------------------------------
    % Package results.
    % ------------------------------------------------------------------------------------------------------------
    PreppedData = struct();
    PreppedData.eeg_epoch_bpf = eeg_epoch_bpf;
    PreppedData.speech_envelope_trials = speech_envelope_trials;
    PreppedData.nfft = nfft;
    PreppedData.num_freqs = num_freqs;
    PreppedData.Freqs = Freqs;
    PreppedData.eeg_Fs = eeg_struct.Fs;
    PreppedData.speech_Fs_downsampled = Speech_Struct.fs; % all trials share the same downsampled rate (== eeg_Fs, by construction)
    PreppedData.Num_channels = Num_channels;
    PreppedData.Num_trials = Num_trials;
    PreppedData.highpass_cutoff = highpass_cutoff;
    PreppedData.lowpass_cutoff = lowpass_cutoff;
    PreppedData.label = eeg_struct.Subj_id;
    PreppedData.Chanlocs = eeg_struct.Chanlocs;
    PreppedData.Condition = eeg_struct.Condition;
    PreppedData.Subj_id = eeg_struct.Subj_id;
    PreppedData.audio_file = eeg_struct.audio_file;
    PreppedData.OnsetLatency_used = OnsetLatency_used;
end
