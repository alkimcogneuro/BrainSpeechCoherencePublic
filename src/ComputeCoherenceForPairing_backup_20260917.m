function [CoherenceResults] = ComputeCoherenceForPairing(PreppedData, pairing, options)
    % ===============================================================================================================
    % COMPUTE CSD / PSD / MSC FOR A GIVEN TRIAL-TO-SPEECH-EPOCH PAIRING
    % ===============================================================================================================
    % Given the pairing-independent preprocessing produced by PreprocessSpeechEEGForCoherence.m, this
    % function computes trial-level CSD/PSD (via the existing, validated CrossSpectralDensity.m --
    % unchanged, called exactly as AnalyzeSpeechEEGCoherence_Dataset.m calls it) and then averages
    % across trials to get CSD_chanmeans / PSD_speech_chanmeans / PSD_eeg_chanmeans / MSC_chanmeans,
    % for whatever trial-to-speech-epoch pairing you specify.
    %
    % "pairing" is a permutation vector of length Num_trials: EEG trial i is analyzed against the
    % speech epoch stored in PreppedData.speech_envelope_trials(:, pairing(i)).
    %   - pairing = 1:Num_trials reproduces the TRUE / observed analysis (trial i's own EEG epoch
    %     paired with trial i's own speech epoch) -- i.e., exactly what AnalyzeSpeechEEGCoherence_Dataset
    %     computes with RandomizeOnsets=false.
    %   - pairing = derangeArray(1:Num_trials) reproduces one draw of the RandomizeOnsets=true control
    %     (every trial's EEG is paired with a DIFFERENT trial's speech epoch). Calling this repeatedly
    %     with fresh derangements is what ShuffleTest_SpeechEEGCoherence.m uses to build a null
    %     distribution.
    %
    % Arguments:
    %   PreppedData: struct returned by PreprocessSpeechEEGForCoherence.m
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
        PreppedData struct
        pairing (1,:) double {mustBePositive, mustBeInteger} = 1:PreppedData.Num_trials
        options.KeepPerTrialValues (1,1) logical = false
    end

    Num_channels = PreppedData.Num_channels;
    Num_trials = PreppedData.Num_trials;
    num_freqs = PreppedData.num_freqs;
    nfft = PreppedData.nfft;
    fs = PreppedData.eeg_Fs;

    if numel(pairing) ~= Num_trials
        error('ComputeCoherenceForPairing:BadPairingLength', ...
            'pairing must have length Num_trials (%d); got length %d.', Num_trials, numel(pairing));
    end
    if numel(unique(pairing)) ~= Num_trials || any(pairing < 1) || any(pairing > Num_trials)
        error('ComputeCoherenceForPairing:BadPairing', ...
            'pairing must be a permutation of 1:%d (each speech epoch used exactly once).', Num_trials);
    end

    csd_vals = zeros(Num_trials, Num_channels, num_freqs);
    psd_speech = zeros(Num_trials, Num_channels, num_freqs);
    psd_eeg = zeros(Num_trials, Num_channels, num_freqs);

    for trial_idx = 1:Num_trials
        speech_envelope = PreppedData.speech_envelope_trials(:, pairing(trial_idx));
        for ch_idx = 1:Num_channels
            this_eeg_epoch_bpf = squeeze(PreppedData.eeg_epoch_bpf(ch_idx, :, trial_idx));
            [csd_vals(trial_idx, ch_idx, :), psd_speech(trial_idx, ch_idx, :), psd_eeg(trial_idx, ch_idx, :)] = ...
                CrossSpectralDensity(speech_envelope, this_eeg_epoch_bpf, nfft, fs);
        end
    end

    CoherenceResults = struct();
    CoherenceResults.CSD_chanmeans = squeeze(mean(csd_vals, 1));
    CoherenceResults.PSD_speech_chanmeans = squeeze(mean(psd_speech, 1));
    CoherenceResults.PSD_eeg_chanmeans = squeeze(mean(psd_eeg, 1));
    CoherenceResults.MSC_chanmeans = abs(CoherenceResults.CSD_chanmeans).^2 ./ ...
        (CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);
    CoherenceResults.Freqs = PreppedData.Freqs;
    CoherenceResults.nfft = nfft;
    CoherenceResults.pairing = pairing;

    if options.KeepPerTrialValues
        CoherenceResults.CSD = csd_vals;
        CoherenceResults.PSD_speech = psd_speech;
        CoherenceResults.PSD_eeg = psd_eeg;
    end
end
