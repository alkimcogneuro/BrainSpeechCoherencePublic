function [CoherenceResults] = ComputeCoherenceForPairing(EEG_Speech_prepped, pairing, options)
    % ===============================================================================================================
    % COMPUTE CSD, PSD, and MSC FOR A GIVEN TRIAL of EEG-SPEECH PAIRING
    % ===============================================================================================================
    % Given the pairing-independent preprocessing produced by PreprocessSpeechEEGForCoherence.m, this
    % function computes trial-level CSD/PSD and then averages across trials to get
    % CSD_chanmeans / PSD_speech_chanmeans / PSD_eeg_chanmeans / MSC_chanmeans, for whatever
    % trial-to-speech-epoch pairing you specify.
    %
    % As of 2026-09-17, each trial's own CSD/PSD estimate comes from WelchCrossSpectralDensity.m
    % (several overlapping windows within the trial, itself averaged before the MSC ratio is ever
    % taken -- see that file's header) rather than a single window spanning the whole trial. This
    % function's own across-TRIAL averaging below is unchanged: CSD/PSD are still averaged across
    % trials before the ratio, exactly as before -- Welch segmentation just adds a second level of
    % spectra-then-ratio averaging (across segments, within each trial) underneath it, using the
    % window length / overlap fraction PreprocessSpeechEEGForCoherence.m already validated for this
    % dataset (EEG_Speech_prepped.WelchWindowSamples / .WelchOverlapFraction).
    %
    % "pairing" is a permutation vector of length Num_trials: EEG trial i is analyzed against the
    % speech epoch stored in EEG_Speech_prepped.speech_envelope_trials(:, pairing(i)).
    %   - pairing = 1:Num_trials reproduces the TRUE / observed analysis (trial i's own EEG epoch
    %     paired with trial i's own speech epoch) -- i.e., exactly what AnalyzeSpeechEEGCoherence_Dataset
    %     computes with RandomizeOnsets=false.
    %   - pairing = derangeArray(1:Num_trials) reproduces one draw of the RandomizeOnsets=true control
    %     (every trial's EEG is paired with a DIFFERENT trial's speech epoch). Calling this repeatedly
    %     with fresh derangements is what ShuffleTest_SpeechEEGCoherence.m uses to build a null
    %     distribution.
    %
    % Arguments:
    %   EEG_Speech_prepped: struct returned by PreprocessSpeechEEGForCoherence.m
    %   pairing:     (1 x Num_trials) permutation of 1:Num_trials. Defaults to the true/observed pairing.
    %   options.KeepPerTrialValues: if true, also return the raw per-trial CSD/PSD_speech/PSD_eeg
    %     arrays ([Num_trials x Num_channels x num_freqs]), matching the .CSD/.PSD_speech/.PSD_eeg
    %     fields of AnalyzeSpeechEEGCoherence_Dataset's CoherenceResults. Default false, since a
    %     shuffle test calls this function hundreds/thousands of times and only needs the
    %     trial-averaged MSC_chanmeans each time -- keeping the per-trial arrays around for every
    %     shuffle would be a lot of wasted memory.
    %
    % Returns CoherenceResults with fields CSD_chanmeans, PSD_speech_chanmeans, PSD_eeg_chanmeans,
    % MSC_chanmeans (all [Num_channels x num_freqs]), Freqs, nfft, pairing (echoed back for
    % bookkeeping), and -- if requested -- CSD / PSD_speech / PSD_eeg per-trial arrays.
    % ===============================================================================================================
    arguments
        EEG_Speech_prepped struct
        pairing (1,:) double {mustBePositive, mustBeInteger} = 1:EEG_Speech_prepped.Num_trials
        options.KeepPerTrialValues (1,1) logical = false
    end

    Num_channels = EEG_Speech_prepped.Num_channels;
    Num_trials = EEG_Speech_prepped.Num_trials;
    num_freqs = EEG_Speech_prepped.num_freqs;
    nfft = EEG_Speech_prepped.nfft;     % Welch window length in samples (see PreprocessSpeechEEGForCoherence.m)
    fs = EEG_Speech_prepped.eeg_Fs;     
    overlap_fraction = EEG_Speech_prepped.WelchOverlapFraction;

    if numel(pairing) ~= Num_trials
        error('ComputeCoherenceForPairing:BadPairingLength', ...
            'pairing must have length Num_trials (%d); got length %d.', Num_trials, numel(pairing));
    end
    if numel(unique(pairing)) ~= Num_trials || any(pairing < 1) || any(pairing > Num_trials)
        error('ComputeCoherenceForPairing:BadPairing', ...
            'pairing must be a permutation of 1:%d (each speech epoch used exactly once).', Num_trials);
    end

    % initialize arrays to hold per-trial CSD and PSD values
    csd_vals = zeros(Num_trials, Num_channels, num_freqs);
    psd_speech = zeros(Num_trials, Num_channels, num_freqs);
    psd_eeg = zeros(Num_trials, Num_channels, num_freqs);

    % loop through trials and channels to compute CSD and PSD for each trial-channel pair
    % use WelchCrossSpectralDensity.m.
    for trial_idx = 1:Num_trials
        speech_envelope = EEG_Speech_prepped.speech_envelope_trials(:, pairing(trial_idx));
        for ch_idx = 1:Num_channels
            eeg_epoch_1chan = squeeze(EEG_Speech_prepped.eeg_epoch_bpf(ch_idx, :, trial_idx));
            [csd_vals(trial_idx, ch_idx, :), psd_speech(trial_idx, ch_idx, :), psd_eeg(trial_idx, ch_idx, :)] = ...
                WelchCrossSpectralDensity(speech_envelope, eeg_epoch_1chan, fs, nfft, overlap_fraction);
        end
    end

    CoherenceResults = struct();
    CoherenceResults.CSD_chanmeans = squeeze(mean(csd_vals, 1));
    CoherenceResults.PSD_speech_chanmeans = squeeze(mean(psd_speech, 1));
    CoherenceResults.PSD_eeg_chanmeans = squeeze(mean(psd_eeg, 1));
    CoherenceResults.MSC_chanmeans = abs(CoherenceResults.CSD_chanmeans).^2 ./ ...
        (CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);
    CoherenceResults.Freqs = EEG_Speech_prepped.Freqs;
    CoherenceResults.nfft = nfft;
    CoherenceResults.pairing = pairing;

    if options.KeepPerTrialValues
        CoherenceResults.CSD = csd_vals;
        CoherenceResults.PSD_speech = psd_speech;
        CoherenceResults.PSD_eeg = psd_eeg;
    end
end
