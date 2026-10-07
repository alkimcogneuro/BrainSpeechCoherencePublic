function [Speech_Struct] = preprocess_speech_epoch(speech_epoch, fs_speech_original, highpass_cutoff, lowpass_cutoff, num_samples_reduced)    
    %-------------------------------------------------------------------------------------------%
    % Prepare an epoch of speech for entrainment anlysis. 
    % Extract the amplitude envelope
    % filter   
    % features of the returned structure are
    % 
    % .envelope: the band-pass filtered speech envelope, downsampled to match EEG sampling rate
    % .phasevals: the instantaneous phase values of the speech envelope
    % .phasevals_rev: the instantaneous phase values of the reversed speech envelope
    % .fs: the new sampling rate of the speech envelope (after downsampling)
    % .fsorig: the original sampling rate of the speech signal
    % .epoch_length: the length of the processed speech envelope vectors.

    % should add optional arguments.
    %-------------------------------------------------------------------------------------------%
    num_samples_speech = length(speech_epoch);
    speech_analytic_signal = hilbert(speech_epoch);                             % compute the analytic signal
    speech_envelope = abs(speech_analytic_signal);                              % extract the amplitude envelope
    %% skipping normalization, bcs we already normalized the raw speech signal before extracting the envelope.  
    %  Normalizing the envelope here would be redundant and could even cause problems if the envelope has very low variance 
    % (e.g., if the speech signal is very quiet or has a lot of silence, which can happen in some epochs).
    %%^ speech_envelope = (speech_envelope - mean(speech_envelope)) / std(speech_envelope);     % normalize
    %% figure; plot(speech_envelope);

    if (any(~isfinite(speech_envelope)))
        % this could happen if the speech vector is all zeros, which has occured before, due to an error in upstream processing.
        error('ERROR . \nSpeech envelope contains non-finite values\n');
    end
    % ==========================================================================================
    % We'll filter the amplitude envelope, to focus on the frequency range of interest for speech processing.
    % The filter parameters will vary, however.
    % For phase coherence, we'll filter between in a band of ~4-8 Hz, 
    % to focus on frequencies that are thought to be relevant for brain-speech entrainment.
    % For Cross Spectral Density, we'll filter in a wider band ~2-35 Hz, to capture the full range of frequencies 
    % that contribute to the speech envelope.  We'll then narrow our focus to the relevant frequency band later,
    % after computing the CSD spectrum.
    % ==========================================================================================
    speech_envelope_bpf = band_pass_filt(speech_envelope, fs_speech_original, highpass_cutoff, lowpass_cutoff);  % filter
    % Downsample from speech sampling rate (e.g,. 44.1 kHz) to the EEG sample rate (e.g., 256 Hz).
    % So that the speech and EEG vectors are the same length            
    %%     I am going to downsample to match the EEG before extracting phase values.
    %%     in the past, I extracted phase first and then downsampled the phase vector.
    %%     I think that was unnecessary and may even have caused problems.
    %% resample() will downsample the speech envelope to match the EEG sampling frequency, 
    % while also applying an anti-aliasing filter.
    speech_envelope_bpf_downsampled = resample(speech_envelope_bpf, num_samples_reduced, num_samples_speech);     % Downsample to match EEG sampling frequency
    Speech_Struct = struct;       % Initialize result variable
    Speech_Struct.envelope = speech_envelope_bpf_downsampled;    % we return the downsampled, band-pass filtered speech envelope, which is what we will use for the entrainment analysis.
    Speech_Struct.phasevals = unwrap(angle(hilbert(speech_envelope_bpf_downsampled)));          % Extract phase from the (downsampled) envelope
    Speech_Struct.phasevals_rev = unwrap(angle(hilbert(flip(speech_envelope_bpf_downsampled))));% Extract phase from the reversed speech envelope 
    % this seems wrong.  can we check it?
    Speech_Struct.fs = fs_speech_original * (num_samples_reduced / num_samples_speech);  % Calculate the new sampling rate for the speech.
    Speech_Struct.fsorig = fs_speech_original;        
    Speech_Struct.epoch_length = length(speech_envelope_bpf_downsampled);       % the length of the processed speech envelope vectors.  this should stay the same....        
end 
