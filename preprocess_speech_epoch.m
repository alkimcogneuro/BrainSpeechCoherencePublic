function [Speech_Struct] = preprocess_speech_epoch(speech_epoch, fs, highpass_cutoff, lowpass_cutoff, num_samples_reduced)    
    %-------------------------------------------------------------------------------------------%
    % Prepare an epoch of speech for entrainment anlysis. 
    % Extract the amplitude envelope
    % filter   
    % features of the returned structure are
    % 
    % envelope: the band-pass filtered speech envelope, downsampled to match EEG sampling rate
    % phasevals: the instantaneous phase values of the speech envelope
    % phasevals_rev: the instantaneous phase values of the reversed speech envelope
    % fs: the new sampling rate of the speech envelope (after downsampling)
    % fsorig: the original sampling rate of the speech signal
    % epoch_length: the length of the processed speech envelope vectors.

    % should add optional arguments.
    %-------------------------------------------------------------------------------------------%
    num_samples_speech = length(speech_epoch);
    speech_analytic_signal = hilbert(speech_epoch);                             % compute the analytic signal
    speech_envelope = abs(speech_analytic_signal);                              % extract the amplitude envelope
    %% skipping normalization, bcs we already normalized the raw speech signal before extracting the envelope.  
    % Normalizing the envelope again would just be redundant and could even cause problems if the envelope has very low variance (e.g., if the speech signal is very quiet or has a lot of silence, which can happen in some epochs).
    %%^ speech_envelope = (speech_envelope - mean(speech_envelope)) / std(speech_envelope);     % normalize
    %% figure; plot(speech_envelope);

    %% plot(speech_envelope);                 % 
    if (any(~isfinite(speech_envelope)))
        % this could happen if the speech vector is all zeros, which has occured before, due to an error in upstream processing.
        error('ERROR . \nSpeech envelope contains non-finite values\n');
    end
    % We'll always filter the amplitude envelope.
    % The filter parameters will vary, however.
    % For phase coherence, we'll filter between 4 and 8 Hz, 
    % to capture the theta band frequencies that are most relevant for speech processing.
    % For Cross Spectral Density, we'll filter between 1 and 35 Hz, to capture the full range of frequencies 
    % that contribute to the speech envelope. After CSD is computed, we will extract the relevant frequency band 
    % from the CSD spectrum for our analysis.
    speech_envelope_bpf = band_pass_filt(speech_envelope, fs, highpass_cutoff, lowpass_cutoff);  % filter
    speech_envelope_bpf_downsampled = resample(speech_envelope_bpf, num_samples_reduced, num_samples_speech);     % Downsample to match EEG sampling frequency
    % Downsample from speech sampling rate (e.g,. 44.1 kHz) to the EEG sample rate (e.g., 256 Hz).
    % So that the speech and EEG vectors are the same length            
    %%     I am going to downsample to match the EEG before extracting phase values.
    %%     in the past, I extracted phase first and then downsampled the phase vector.
    %%     i think that was unnecessary and may even have caused problems.
    Speech_Struct = struct;       % Initialize result variable
    Speech_Struct.envelope = speech_envelope_bpf_downsampled;
    Speech_Struct.phasevals = unwrap(angle(hilbert(speech_envelope_bpf_downsampled)));          % Extract phase from the envelope
    Speech_Struct.phasevals_rev = unwrap(angle(hilbert(flip(speech_envelope_bpf_downsampled))));% Extract phase from the reversed speech envelope 
    Speech_Struct.fs = fs * (num_samples_reduced / num_samples_speech);  % Calculate the new sampling rate for the speech.
    Speech_Struct.fsorig = fs;        
    Speech_Struct.epoch_length = length(speech_envelope_bpf_downsampled);       % the length of the processed speech envelope vectors.  this should stay the same....        
end 
