function [CoherenceResults] = Apply2Dataset_InstaPhaseCoherence(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % This function applies Instantaneous Phase Coherence analysis to one EEG dataset. 
    %
    % Arguments:
    %   eeg_struct:     structure containing EEG data and metadata
    %   speech_rawdata: structure containing raw speech data and metadata
    %   highpass_cutoff: high-pass filter cutoff frequency (Hz)
    %   lowpass_cutoff:  low-pass filter cutoff frequency (Hz)
    %
    % example call with optional arguments provided
    % Apply2Dataset_InstantPhaseCoher(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, StartTimeOffset=0.2, EpochDuration=1.8)
    % Returns:
    %   CoherenceResults: structure containing results
    % 
    % NOTE:  the highpass_cutoff and lowpas_cutoff parameters here mean something different from the same parameters CSD functions
    %  
    
    % ===============================================================================================================
    arguments
        eeg_struct
        speech_rawdata
        highpass_cutoff (1,1) double {mustBeReal} = 2
        lowpass_cutoff  (1,1) double {mustBeReal} = 35
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.  
                                                                % If NaN, use the original EEG onset as the start of the analysis epoch.
        options.EpochDuration (1,1) double {mustBeReal} = NaN   % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
        options.RandomizeOnsets (1,1) logical = false  % optional boolean argument to indicate whether to randomize the speech onset latencies for the analysis, as a control.  
    end
    fprintf(' -- Running Instantaneous Phase Coherence analysis for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Filter parameters: HighPass = %d Hz, LowPass = %d Hz\n', highpass_cutoff, lowpass_cutoff);
    %% fprintf(' --Optional parameters: StartTimeOffset = %.2f seconds, EpochDuration = %.2f seconds\n', options.StartTimeOffset, options.EpochDuration);
    fprintf(' -- EEG sampling rate: %d Hz, Speech sampling rate: %d Hz\n', eeg_struct.Fs, speech_rawdata.Fs);
    fprintf(' -- Number of EEG trials: %d, Number of EEG channels: %d\n', eeg_struct.Num_trials, eeg_struct.Num_channels);
    fprintf(' -- Length of EEG epochs (raw data): %d samples (%.2f seconds)\n', size(eeg_struct.Data, 2), size(eeg_struct.Data, 2) / eeg_struct.Fs);    
    fprintf(' -- RandomizeOnsets: %d\n', options.RandomizeOnsets);
       

    
    %--------------------------------------------------------------------------------------------------------------------
    % We will check the validity of the optional parameters for epoching, and if they are not valid, throw an error to alert the user to adjust them:
    % 1. If either StartTimeOffset or EpochDuration is provided, then both must be provided, 
    %   because we need to know both the offset and the duration in order to determine the new speech epoch start and end times,
    %   which will be aligned with the new EEG epoch.  If we adjust one without the other, then the speech and EEG epochs will be misaligned in time.
    % 2. The StartTimeOffset must be non-negative.
    % 3. The EpochDuration must be less than or equal to the length of the EEG epochs, 
    %    because if the epoch duration is longer than the length of the EEG epochs, then we will be trying to analyze data that doesn't exist, which will lead to errors in the analysis.
    % 4. The StartTimeOffset must be less than or equal to the EpochDuration, because if the offset is greater than the epoch duration, 
    %    then we'd be asking to offset the speech onset to a time that is after the end of the EEG epoch, which would lead to misalignment between the speech and EEG epochs.
    % 5. If the optional parameters are not provided (i.e., they are NaN), then we will use the full length of the EEG epochs as the analysis epoch duration, and we will not apply any offset to the speech onset latencies, so that the speech and EEG epochs will be aligned in time based on the original EEG onset.
    % 6. We should also check that the StartTimeOffset and EpochDuration are not too long for the length of the EEG epochs, and if they are, throw an error to alert the user to adjust them or use NaN to use the full length of the EEG epochs.
    %--------------------------------------------------------------------------------------------------------------------
    if ~isnan(options.StartTimeOffset) || ~isnan(options.EpochDuration)
        fprintf(' -- Optional parameters for epoching provided:\n');    
        % if either StartTimeOffset or EpochDuration is provided, then both must be provided, because we need to know both the offset and the duration in order to determine the new speech epoch start and end times, which will be aligned with the new EEG epoch.  If we adjust one without the other, then the speech and EEG epochs will be misaligned in time.
        if ~isnan(options.StartTimeOffset) && isnan(options.EpochDuration)
            error('If StartTimeOffset is provided, EpochDuration must also be provided, because we need to know the duration of the analysis epoch in order to determine the new speech epoch start and end times, which will be aligned with the new EEG epoch.');
        end
        if ~isnan(options.EpochDuration) && isnan(options.StartTimeOffset)
            error('If EpochDuration is provided, StartTimeOffset must also be provided, because if we adjust the epoch duration without adjusting the start time, then the speech and EEG epochs will be misaligned in time.');
        end
        % The StartTimeOffset must be non-negative.
        if ~isnan(options.StartTimeOffset) && (options.StartTimeOffset < 0)
            error('The specified StartTimeOffset of %.2f seconds is negative. Please adjust the StartTimeOffset to be non-negative, or use NaN to use the original EEG onset as the start of the analysis epoch.', options.StartTimeOffset);
        end
        % The EpochDuration must be less than or equal to the length of the EEG epochs, because if the epoch duration is longer than the length of the EEG epochs, then we will be trying to analyze data that doesn't exist, which will lead to errors in the analysis.
        if ~isnan(options.EpochDuration) && (options.EpochDuration > (size(eeg_struct.Data, 2) / eeg_struct.Fs))
            error('The specified EpochDuration of %.2f seconds is too long for the length of the EEG epochs, which is %.2f seconds. Please adjust the EpochDuration to be less than or equal to the length of the EEG epochs, or use NaN to use the full length of the EEG epochs.', options.EpochDuration, size(eeg_struct.Data, 2) / eeg_struct.Fs);
        end
        % The StartTimeOffset must be less than or equal to the EpochDuration, because if the offset is greater than the epoch duration, then we'd be asking to offset the speech onset to a time that is after the end of the EEG epoch, which would lead to misalignment between the speech and EEG epochs.
        if ~isnan(options.StartTimeOffset) && ~isnan(options.EpochDuration) && (options.StartTimeOffset > options.EpochDuration)
             error('The specified StartTimeOffset of %.2f seconds is greater than the specified EpochDuration of %.2f seconds. Please adjust the StartTimeOffset to be less than or equal to the EpochDuration, because if the offset is greater than the epoch duration, it offsets the speech onset to a time that is after the end of the EEG epoch, which would lead to misalignment between the speech and EEG epochs.', options.StartTimeOffset, options.EpochDuration);
        end
        % If the user provided a StartTimeOffset and EpochDuration, then we will 
        % adjust the EEG data structure to only include the samples in the specified range relative to the original EEG onset.  
        % This will shift the analysis epoch for both the EEG and speech data by the specified offset and duration.  
        % If the user did not provide these optional arguments, 
        % then we will use the full length of the EEG epochs as the analysis epoch duration. 
        if (~isnan(options.StartTimeOffset) & ~isnan(options.EpochDuration))
            % User has provided a StartTimeOffset and EpochDuration, 
            % so we will adjust the EEG data structure to only include the samples in the specified range relative to the original EEG onset.
            % if both the optional parameters are provided correctly, then print them for display
            fprintf('    - StartTimeOffset: %.2f seconds, which corresponds to %d samples at the EEG sampling rate\n', options.StartTimeOffset, round(options.StartTimeOffset * eeg_struct.Fs));
            fprintf('    - EpochDuration: %.2f seconds, which corresponds to %d samples at the EEG sampling rate\n', options.EpochDuration, round(options.EpochDuration * eeg_struct.Fs));    
            % Round and convert the start time offset from seconds to samples, 
            % so that we can apply it to the speech onset latencies in the EEG data structure.
            start_offset_samples = round(options.StartTimeOffset * eeg_struct.Fs); % the number of samples corresponding to the specified StartTimeOffset, at the EEG sampling rate.  
            epoch_duration_samples = round(options.EpochDuration * eeg_struct.Fs); % the number of samples corresponding to the specified EpochDuration, at the EEG sampling rate. 
            sample_range_restricted = start_offset_samples:(start_offset_samples + epoch_duration_samples - 1);  % The range of samples to use for the analysis epoch, relative to the original EEG onset.  This will be applied to the speech onset latencies to determine the new speech epoch start and end times.
            fprintf(' -- Applying StartTimeOffset of %.2f seconds and EpochDuration of %.2f seconds to EEG data structure.\n', options.StartTimeOffset, options.EpochDuration);
            % Restrict the EEG data to the specified sample range for the analysis epoch.
            eeg_struct.Data = eeg_struct.Data(:, sample_range_restricted, :);
            % shift the speech onset latency by the specified offset, so that the new speech epoch will be aligned with the new EEG epoch. 
            % this will work for all trials, because the speech epoch will be defined relative to the new EEG onset for each trial.
            % eeg_struct.OnsetLatency is a vector of length Num_trials, with the original speech onset latency for each trial.
            % we will add the StartTimeOffset (in seconds) to each of these latencies, to get the new speech onset latency for each trial, which will be aligned with the new EEG epoch. 
            % Note:  if the StartTimeOffset is present, we will adjust the speech epoch duration to match the new EEG epoch duration, 
            % so that the speech and EEG epochs will be aligned in time.
            eeg_struct.OnsetLatency = eeg_struct.OnsetLatency + options.StartTimeOffset;
        end
    else
        fprintf(' -- No optional parameters for epoching provided; using full length of EEG epochs for analysis.\n');
    end

    % The data matrix should be [num_channels x num_samples x num_trials]
    Num_channels = size(eeg_struct.Data, 1);  
    Num_samples = size(eeg_struct.Data, 2);  % the number of samples in the EEG epochs after any optional adjustments to the epoch duration.  This will be used to determine the duration of the speech epochs and to set the nfft parameter for the CSD calculation.
    Num_trials = size(eeg_struct.Data, 3);

    % Z-score normalization of EEG data across all samples and trials, separately for each channel.
    % Reshape to (numchannels x (numsamples*numtrials)) in order to copmute mean and std across all samples and trials for each channel, 
    % then reshape back to original dimensions after normalization.

    % ****
    %    For Phase Coherence, we don't need to normalize the EEG or the speech data, 
    % because we are only interested in the phase information, 
    % and the phase is not affected by the amplitude of the signal.
    % I think we should get rid of these normalization steps for the Instantaneous Phase Coherence analysis, 
    % because they are not necessary and could potentially introduce artifacts into the phase information.

    eegdata_2d = reshape(eeg_struct.Data, Num_channels, []);
    mu    = mean(eegdata_2d, 2);        % 64x1 mean per channel
    sigma = std(eegdata_2d, 0, 2);      % 64x1 std per channel
    eeg_data_norm = (eegdata_2d - mu) ./ sigma;    % Subtract mean and divide by std  
    eeg_struct.Data = reshape(eeg_data_norm, Num_channels, Num_samples, Num_trials);  % reshape back to original dimensions
    % Z-score normalization of speech data.
    % We normalize the speech across the entire speech signal
    % because the speech is not segmented into trials in the same way as the EEG data, and we want to preserve the relative amplitude differences across the entire speech signal,
    speech_amplitudes = speech_rawdata.Amplitudes;
    speech_amplitudes_norm = (speech_amplitudes - mean(speech_amplitudes)) / std(speech_amplitudes);
    speech_rawdata.Amplitudes = speech_amplitudes_norm;  % update the speech raw data structure with the normalized amplitudes.
    
    % nfft is the number of points to use in the FFT calculation for the CSD estimation.
    % we want to set nfft equal to the number of samples in the EEG epochs, so that we are computing the FFT on the full length of the analysis epochs.
    % note that if nfft is a power of 2 for computational efficiency; 
    % this will slightly reduce the frequency resolution but is generally a good idea for FFT calculations.
    % for example, if the EEG epochs are 1500 samples long, then nfft will be set to 1024, which is the largest power of 2 less than or equal to 1500. 
    % This means that the FFT will be computed on the first 1024 samples of the EEG and speech epochs, and the remaining samples will be ignored for the CSD calculation. 
    % nfft = 2^floor(log2(Num_samples));
    nfft = Num_samples;
    % find the largest even number of samples that is less than or equal to Num_samples, to ensure that nfft is even, 
    % which is important for interpreting the frequency bins of the FFT output.
    %if mod(nfft, 2) ~= 0
    %    nfft = nfft - 1;  % make nfft even if it is odd
    %end
    
    fprintf(' --Number of samples in EEG epochs: %d samples (%.2f seconds)\n', Num_samples, (Num_samples / eeg_struct.Fs));
    fprintf(' --Number of FFT points (nfft) for CSD calculation: %d\n', nfft);  
    
    CoherenceResults = struct();  % Initialize
    CoherenceResults.analysis_type = 'Instantaneous_Phase_Coherence';     
    CoherenceResults.highpass_cutoff = highpass_cutoff;    
    CoherenceResults.lowpass_cutoff = lowpass_cutoff;      
    CoherenceResults.eeg_Fs = eeg_struct.Fs;                               
    CoherenceResults.speech_Fs = speech_rawdata.Fs;                       
    CoherenceResults.label = eeg_struct.Subj_id;
    CoherenceResults.Chanlocs = eeg_struct.Chanlocs;  
    % Save additional information about the EEG data and analysis parameters in the main results structure, for reference.
    CoherenceResults.Block = eeg_struct.Block;
    CoherenceResults.exposure_type = eeg_struct.exposure_type;
    CoherenceResults.Subj_id = eeg_struct.Subj_id;
    CoherenceResults.audio_file = eeg_struct.audio_file;
    CoherenceResults.nfft = nfft;
    CoherenceResults.Num_channels = Num_channels;
    % Initialize results matrix: NumTrials x NumChannels
    cross_spectral_density_vals = cell(eeg_struct.Num_trials, eeg_struct.Num_channels);
    speech_epoch_duration = size(eeg_struct.Data, 2) / eeg_struct.Fs;  % speech epoch duration in seconds, calculated from the number of samples in the EEG epochs and the EEG sampling rate.
    
    % ------------------------------------------------------------------------------------------------------- %   
    % Control Analysis.
    % If the RandomizeOnsets option is true, then we will add a random offset to the speech onset latency for each trial, 
    % to shift the speech epoch to a different segment of the speech signal that does not correspond to the EEG data for that trial, 
    % as a control analysis.
    % ------------------------------------------------------------------------------------------------------- %   
    if options.RandomizeOnsets
        fprintf(' -- Randomized onset mode ON: adding random 3-10s offset to speech onset latencies for each trial, to shift speech epochs to different segments of the speech signal that do not correspond to the EEG data for those trials.\n'); 
        for eeg_trial_idx = 1:eeg_struct.Num_trials
            % For each trial, we will randomly select a different trial index 
            % that is at least 2 trials away from the current trial index (eeg_trial_idx).
            rnd_mismatch_eeg_trial_idx = randi(eeg_struct.Num_trials);
            while abs(rnd_mismatch_eeg_trial_idx - eeg_trial_idx) < 2  
                % loop until we get a trial index that satisfies our critrion.
                rnd_mismatch_eeg_trial_idx = randi(eeg_struct.Num_trials);
            end
            eeg_struct.rand_mismatch_indices(eeg_trial_idx) = rnd_mismatch_eeg_trial_idx;  % save the randomly selected mismatching trial index in a new field in the EEG data structure, for reference.
            % Store the onset latency for this randomly selected trial in a new field in the EEG data structure called RandOnsetLatency, 
            % which will hold the randomized speech onset latencies for each trial.
            % set the speech onset latency for this trial to the speech onset latency of the randomly selected trial, to shift the speech epoch to a different segment of the speech signal that does not correspond to the EEG data for this trial. 
            RandOnsetLatency(eeg_trial_idx) = eeg_struct.OnsetLatency(rnd_mismatch_eeg_trial_idx);  
        end

        % Save the original speech onset latencies, for record-keeping.  
        OnsetLatency_Original = eeg_struct.OnsetLatency; 
        % Overwrite the speech onset latencies in the EEG data structure with the randomized ones, 
        % so that the rest of the analysis will use the randomized onset latencies.         
        eeg_struct.OnsetLatency = RandOnsetLatency; 
        % Add a field to the main results structure to indicate that we used randomized onsets for this analysis, 
        % and save the original speech onset latencies in the results structure as well, so that we have a record of the original and randomized speech onset latencies for each trial in the results of this analysis, which will be important for interpreting the results of the control analysis and comparing them to the main analysis with the original speech onset latencies.
        CoherenceResults.RandomizedOnsets = true;
        CoherenceResults.OnsetLatency_Original = OnsetLatency_Original; % save the original speech onset latencies in the results structure for reference.
        CoherenceResults.OnsetLatency_Randomized = RandOnsetLatency;    % save the randomized speech onset latencies in the results structure for reference.
    end
    % ----------------------------------------------------------------------------------------------------------------- %
    % Loop through trials to extract speech epochs and then compute phase coherence analysis
    % with EEG data for each trial and channel.
    % ----------------------------------------------------------------------------------------------------------------- %
    for eeg_trial_idx = 1:eeg_struct.Num_trials
        speech_onset_latency = eeg_struct.OnsetLatency(eeg_trial_idx);          % speech onset latency for this trial, in seconds.
        speech_offset_latency = speech_onset_latency + speech_epoch_duration;   % speech offset latency for this trial, in seconds
        speech_onset_idx = uint64(speech_onset_latency * speech_rawdata.Fs);    % Convert seconds to samples 
        speech_offset_idx = uint64(speech_offset_latency * speech_rawdata.Fs);  % Convert seconds to samples
        % Ensure indices are within bounds
        speech_onset_idx = max(1, min(speech_onset_idx, length(speech_rawdata.Amplitudes)));
        speech_offset_idx = max(1, min(speech_offset_idx, length(speech_rawdata.Amplitudes)));
        speech_epoch = speech_rawdata.Amplitudes(speech_onset_idx:speech_offset_idx);  % extract the speech epoch.
        % Prepare the speech for coherence analysis
        % Extract amplitude envelope, bandpass filter, downsample, save phase and magnitude information. 
        % I was previously assuming that eeg_struct would have a field called Num_samples, but it doesn't,
        % so I'm just calculating the number of samples from the size of the data matrix.
        num_samples = size(eeg_struct.Data, 2);  % Number of samples in the EEG epoch; we will downsample the speech to match this number of samples.
        fprintf('run preprocessing for speech epoch: trial %d, speech onset latency = %.2f seconds, speech offset latency = %.2f seconds, number of samples in speech epoch = %d\n', eeg_trial_idx, speech_onset_latency, speech_offset_latency, length(speech_epoch));
        % Speech_struct will contain the preprocessed speech information that we will use for the coherence analysis, 
        % including the amplitude envelope, phase values, and magnitude values for each frequency bin.
         
        [Speech_Struct] = preprocess_speech_epoch(speech_epoch, speech_rawdata.Fs, ...
                                                  highpass_cutoff, lowpass_cutoff, num_samples);
        % Loop through EEG channels
        for ch_idx = 1:eeg_struct.Num_channels
            eeg_epoch = squeeze(eeg_struct.Data(ch_idx, :, eeg_trial_idx));   % one epoch of EEG data for this trial and channel, as a vector of length Num_samples.
            eeg_epoch = eeg_epoch(:);  % convert to column vector if it is not already a column vector.
            eeg_epoch_bpf = band_pass_filt(eeg_epoch, eeg_struct.Fs, highpass_cutoff, lowpass_cutoff);  % bandpass filter EEG epoch
            % Compute instantaneous phase of the EEG signal at each time point in the epoch, 
            % using the Hilbert transform to get the analytic signal, 
            % and then taking the angle of the analytic signal to get the instantaneous phase values.
            % unwrap() is used to remove phase discontinuities and provide a continuous phase signal over time, 
            % which is important for accurately calculating the phase differences and coherence values.
            eeg_phasevals = unwrap(angle(hilbert(eeg_epoch_bpf)));
            % Phase difference: EEG vs. speech
            % calculated at each time point in the epoch as the difference 
            % between the instantaneous phase of the EEG signal 
            % and the instantaneous phase of the speech envelope at each time point.
            phase_diffs = eeg_phasevals - Speech_Struct.phasevals;
            % ------------------------------------------------------------------------
            % Phase coherence.
            % We calculate the phase coherence values for each trial and channel 
            % as the magnitude of the mean vector of the phase differences across all time points in the epoch. 
            % This average gives us a measure of how consistent the phase relationship is 
            % between the EEG signal and the speech envelope across time within that trial and channel.
            % Store the phase coherence values for this trial and channel in the results matrix [num_trials x num_channels]. 
            % ------------------------------------------------------------------------
            phase_coherence_vals(eeg_trial_idx, ch_idx) = abs(mean(exp(1i * phase_diffs)));            
        end
    end

    % Save trial-level values in the main results structure.
    CoherenceResults.phase_coherence = phase_coherence_vals;   
    % Average across trials per channel.  
    CoherenceResults.phase_coherence_chanmeans = mean(phase_coherence_vals, 1);
    % Maximum coherence and corresponding channel
    [CoherenceResults.max_val, CoherenceResults.max_chan]  = max(CoherenceResults.phase_coherence_chanmeans);
end

