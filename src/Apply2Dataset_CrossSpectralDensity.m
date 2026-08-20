function [CoherenceResults] = Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % This function applies Cross Spectral Density analysis to one EEG dataset.
    % 
    % Arguments:
    %   eeg_struct:     structure containing EEG data and metadata, including pointers to speech file and onset latencies of critical speech events. 
    %   speech_rawdata: structure containing raw speech data (full contents of one speech file) and metadata
    %   highpass_cutoff: high-pass filter cutoff frequency (Hz)
    %   lowpass_cutoff:  low-pass filter cutoff frequency (Hz)
    %   optional, StartTimeOffset:  If user wants the analysis epoch to begin after the onsets recorded in the EEG Data atructure
    %                               then we provide that information in this optional argument, in seconds. 
    %   
    %   optional, EpochDuration:    the duration, in seconds, of the optinal 
    % 
    % Generally, this will be the EEG data from one experimental subject, in one condition. 
    % example call with optional arguments provided
    % 
    % Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, StartTimeOffset=0.2, EpochDuration=1.8)
    % Returns:
    %   CoherenceResults: structure containing results
    % 
    % Returns a structure CoherenceResults, which contains the following fields
    %  --  .CSD_chanmeans: a [num_channels X num_frequencies] matrix containing the average CSD values across trials, for each channel and frequency bin.
    %  --  .PSD_speech_chanmeans: a [num_channels X num_frequencies] matrix containing the average power spectral density values for the speech signal across trials, for each channel and frequency bin.
    %  --  .PSD_eeg_chanmeans: a [num_channels X num_frequencies] matrix containing the average power spectral density values for the EEG signal across trials, for each channel and frequency bin.
    %  --  .MSC_chanmeans: a [num_channels X num_frequencies] matrix containing the average magnitude squared coherence values across trials, for each channel and frequency bin.
    %                      see below for details on how we calculate MSC from the CSD and PSD values.
    %  --  .analysis_type: a string indicating the type of analysis performed (e.g., 'cross_spectral_density'; different analysis types are returned by other functions in this suite). 
    %  --  .highpass_cutoff: the high-pass filter cutoff frequency used in the analysis, in Hz.
    %  --  .lowpass_cutoff: the low-pass filter cutoff frequency used in the analysis, in Hz.
    %  --  .eeg_Fs: the sampling rate of the EEG data, in Hz.
    %  --  .speech_Fs: the sampling rate of the speech data, in Hz.
    %  --  .label: a string label for the dataset (e.g., subject ID), taken from the eeg_struct.Subj_id field.
    %  --  .Chanlocs: the channel location information from the EEG data structure, which can be used for topographical plotting of the results.
    %  --  .Condition: the experimental condition for this dataset, taken from the eeg_struct.Condition field.
    %  --  .Subj_id: the subject ID for this dataset, taken from the eeg_struct.Subj_id field.
    %  --  .audio_file: the name of the audio file used in this dataset, taken from the eeg_struct.audio_file field.
    %  --  .nfft: the number of points used in the FFT calculation for the CSD estimation, which is determined based on the length of the EEG epochs and any optional adjustments to the epoch duration.
    %  --  .Num_channels: the number of EEG channels in the dataset, taken from the size of the eeg_struct.Data matrix.
    % 
    % 
    % Example usage:
    %    eeg_struct = load_eeg_data('subject01_conditionA.mat');  %    Load EEG data structure for one subject and condition
    %    speech_rawdata = load_speech_data('subject01_conditionA_speech.wav');  %    Load raw speech data structure for the corresponding speech file
    %    highpass_cutoff = 2;  %    High-pass filter cutoff frequency in Hz
    %    lowpass_cutoff = 35;  %    Low-pass filter cutoff frequency in Hz
    %    options.StartTimeOffset = 0.2;  %    Optional: start the analysis epoch 0.2 seconds after the original EEG onset
    %    options.EpochDuration = 1.8;    %    Optional: set the analysis epoch duration to 1.8 seconds
    %    CoherenceResults = Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options);
    %    The function will return the CoherenceResults structure containing the CSD, PSD, and MSC results, as well as information about the analysis parameters and the EEG data for reference.

    % ===============================================================================================================
    arguments
        eeg_struct                                      % 
        speech_rawdata
        highpass_cutoff (1,1) double {mustBeReal} = 2
        lowpass_cutoff  (1,1) double {mustBeReal} = 35
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.  
                                                                % If NaN, use the original EEG onset as the start of the analysis epoch.
        options.EpochDuration (1,1) double {mustBeReal} = NaN   % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
        options.RandomizeOnsets (1,1) logical = false  % optional boolean argument to indicate whether to randomize the speech onset latencies for the analysis, as a control.  
    end
        
    fprintf(' -- Running Cross Spectral Density analysis for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Speech file is %s\n', eeg_struct.audio_file);
    fprintf(' -- Filter parameters: HighPass = %.2f Hz, LowPass = %.2f Hz\n', highpass_cutoff, lowpass_cutoff);
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

    % Check EEG_struct.OnsetLatency to make sure that all onset latencies are real numbers and not NaNs.
    if any(isnan(eeg_struct.OnsetLatency))
        error('EEG_struct.OnsetLatency contains NaN values. Please check the EEG data structure and ensure that all onset latencies are valid real numbers.');
    end
    

    % The data matrix should be [num_channels x num_samples x num_trials]
    Num_channels = size(eeg_struct.Data, 1);  
    Num_samples_eeg = size(eeg_struct.Data, 2);  % the number of samples in the EEG epochs after any optional adjustments to the epoch duration.  This will be used to determine the duration of the speech epochs and to set the nfft parameter for the CSD calculation.
    Num_trials = size(eeg_struct.Data, 3);
    fprintf(' -- Number of samples in EEG epochs: %d samples (%.2f seconds)\n', Num_samples_eeg, (Num_samples_eeg / eeg_struct.Fs));
    
    % If the EEG data has more than 64 channels, get rid of all channels beyond the first 64, 
    % to match the number of channels in the other datasets and to focus on the standard 64-channel montage for this analysis.
    % NOTE:  we should change the code to be more flexible about the number of channels, 
    % but for now, we will restrict to 64 channels to match the other datasets and the standard montage.
    if Num_channels > 64
        fprintf(' -- EEG data has %d channels, which is more than 64. Restricting analysis to the first 64 channels to match the standard montage and other datasets.\n', Num_channels);
        eeg_struct.Data = eeg_struct.Data(1:64, :, :);  % keep only the first 64 channels
        eeg_struct.Chanlocs = eeg_struct.Chanlocs(1:64);  % keep only the channel location information for the first 64 channels
        eeg_struct.Num_channels = 64;  % update the Num_channels field in the EEG data structure to reflect the restriction to 64 channels.
        Num_channels = 64;  % update the number of channels variable to reflect the restriction to 64 channels.
    end
    % ------------------------------------------------------------------------------------------------------------------- %
    % Z-score normalization of EEG data across all samples and trials, separately for each channel.
    % Reshape to (numchannels x (numsamples*numtrials)) in order to copmute mean and std across all samples and trials for each channel, 
    % then reshape back to original dimensions after normalization.
    % ------------------------------------------------------------------------------------------------------------------- %
    eegdata_2d = reshape(eeg_struct.Data, Num_channels, []);
    mu    = mean(eegdata_2d, 2);        % 64x1 mean per channel
    sigma = std(eegdata_2d, 0, 2);      % 64x1 std per channel
    eeg_data_norm = (eegdata_2d - mu) ./ sigma;    % Subtract mean and divide by std  
    eeg_struct.Data = reshape(eeg_data_norm, Num_channels, Num_samples_eeg, Num_trials);  % reshape back to original dimensions
    % Z-score normalization of speech data.
    % We normalize the speech across the entire speech signal
    % because the speech is not segmented into trials in the same way as the EEG data, and we want to preserve the relative amplitude differences across the entire speech signal,
    speech_amplitudes = speech_rawdata.Amplitudes;
    speech_amplitudes_norm = (speech_amplitudes - mean(speech_amplitudes)) / std(speech_amplitudes);
    speech_rawdata.Amplitudes = speech_amplitudes_norm;  % update the speech raw data structure with the normalized amplitudes.
    speech_rawdata.length_sec = length(speech_rawdata.Amplitudes) / speech_rawdata.Fs;  % length of the speech in seconds 
    fprintf(' -- length of the speech vector in seconds: %.2f', speech_rawdata.length_sec)
    % nfft is the number of points to use in the FFT calculation for the CSD estimation.
    % We will set nfft equal to the number of samples in the EEG epochs, so that we are computing the FFT on the full length of the analysis epochs.
    nfft = Num_samples_eeg;    
    
    fprintf(' -- Number of FFT points (nfft) for CSD calculation: %d\n', nfft);  
    CoherenceResults = struct();  % Initialize
    CoherenceResults.analysis_type = 'cross_spectral_density';     
    CoherenceResults.highpass_cutoff = highpass_cutoff;    
    CoherenceResults.lowpass_cutoff = lowpass_cutoff;      
    CoherenceResults.eeg_Fs = eeg_struct.Fs;                               
    CoherenceResults.speech_Fs = speech_rawdata.Fs;                       
    CoherenceResults.label = eeg_struct.Subj_id;
    CoherenceResults.Chanlocs = eeg_struct.Chanlocs;  
    % Save additional information about the EEG data and analysis parameters in the main results structure, for reference.
    CoherenceResults.Condition = eeg_struct.Condition;
    CoherenceResults.Subj_id = eeg_struct.Subj_id;
    CoherenceResults.audio_file = eeg_struct.audio_file;
    CoherenceResults.nfft = nfft;
    CoherenceResults.Num_channels = Num_channels;
    speech_epoch_duration_seconds = size(eeg_struct.Data, 2) / eeg_struct.Fs;  % speech epoch duration in seconds, calculated from the number of samples in the EEG epochs and the EEG sampling rate.
    speech_epoch_duration_samples = round(speech_epoch_duration_seconds * speech_rawdata.Fs);  % convert speech duration from seconds to samples (at the speech sample rate)
                                                                                            % we round, in case the result is non-integer.
    % ------------------------------------------------------------------------------------------------------- %   
    % CONTROL ANALYSIS.
    % If the RandomizeOnsets option is true, then we will add a random offset to the speech onset latency for each trial, 
    % to shift the speech epoch to a different segment of the speech signal that does not correspond to the EEG data for that trial, 
    % as a control analysis.
    % ------------------------------------------------------------------------------------------------------- %   
    if options.RandomizeOnsets
        fprintf(' -- Randomized onset mode ON: adding 4s offset to speech onset latencies for each trial, to shift speech epochs to different segments of the speech signal that do not correspond to the EEG data for those trials.\n'); 
        RandOnsetLatency = derange(eeg_struct.OnsetLatency);  % scramble the order of the original onset latencies, every trial has a the onset associated with a different trial.
 
        %{
         for eeg_trial_idx = 1:eeg_struct.Num_trials
            % To randomize the speech epoch associated with each EEG trial, we will shift the speech onset latency for each trial 
            % by borrowing the onset latency from a different trial, 
            % We'll collect those new onset latencies in RandOnsetLatency, which will be the same length as the number of trials, 
            % and then we will overwrite the OnsetLatency field in the EEG data structure with these new randomized onset latencies,
            % so that the rest of the analysis will use the randomized onset latencies.           
            if (eeg_trial_idx >= 2) 
                % set the randomized onset latency for this trial to be the same as the original onset latency 
                % for the previous trial, so that the speech epoch for this trial will be aligned with the EEG data from the previous trial, 
                % rather than the current trial.
                RandOnsetLatency(eeg_trial_idx) = eeg_struct.OnsetLatency(eeg_trial_idx-1);              
            else
                % but if this is the first trial, then we don't have a previous trial to borrow the onset latency from, 
                % so we will borrow from a downstream trial instead, to shift it forward by 2 trials.
                RandOnsetLatency(eeg_trial_idx) = eeg_struct.OnsetLatency(eeg_trial_idx+2);  % for the first trial, we don't have a previous trial to borrow the onset latency from, so we will just shift it forward by 4 seconds. 
            end
            %-------------------------------------------------------------------------------%
            % For each trial, store a speech onset latency value that is 4 seconds away from the actual onset.
            % Place it 4 seconds before the actual onset, unless that would be within the first 2 seconds of the recording
            % Otherwise, place it 4 seconds after the actual onset. 
            % current_latency = eeg_struct.OnsetLatency(eeg_trial_idx);
            % if ((current_latency - 4) < 2)
            %     RandOnsetLatency(eeg_trial_idx) = current_latency + 4;   % shift forward by 4 seconds
            % else
            %     RandOnsetLatency(eeg_trial_idx) = current_latency - 4;  % shift backward by 4 seconds.
            % end
            %-------------------------------------------------------------------------------%
        end 
%}

        % Save the original speech onset latencies, for record-keeping.  
        OnsetLatency_Original = eeg_struct.OnsetLatency; 
        % Overwrite the speech onset latencies in the EEG data structure with the randomized ones, 
        % so that the rest of the analysis will use the randomized onset latencies.         
        eeg_struct.OnsetLatency = RandOnsetLatency; 
        % Add a field to the main results structure to indicate that we used randomized onsets for this analysis, 
        % and save the original speech onset latencies in the results structure as well, so that we have a record of the original and randomized speech onset latencies for each trial in the results of this analysis, which will be important for interpreting the results of the control analysis and comparing them to the main analysis with the original speech onset latencies.
        CoherenceResults.RandomizedOnsetsFlag = true;
        CoherenceResults.OnsetLatency_Original = OnsetLatency_Original; % save the original speech onset latencies in the results structure for reference.
        CoherenceResults.OnsetLatency_Randomized = RandOnsetLatency;    % save the randomized speech onset latencies in the results structure for reference.
    else 
        CoherenceResults.RandomizedOnsetsFlag = false;
    end    
    % ----------------------------------------------------------------------------------------------------------------- %
    % Loop through trials to extract speech epochs and then
    % compute CSD analysis with EEG data for each trial and channel.
    % ----------------------------------------------------------------------------------------------------------------- %
    % nfft is the number of bins in the one-sided spectrum that CrossSpectralDensity will return.
    % This must match the parity-aware logic inside CrossSpectralDensity itself: nfft/2 + 1
    % bins for even nfft (DC through Nyquist), or (nfft+1)/2 bins for odd nfft (DC through
    % the highest unique bin; there is no exact Nyquist bin when nfft is odd).
    if mod(nfft, 2) == 0
        num_freqs = nfft/2 + 1;
    else
        num_freqs = (nfft + 1)/2;
    end
    % Initialize, for pre-allocation, these three results structures for every iteration of the loop.
    csd_vals = zeros(Num_trials, Num_channels, num_freqs);      % csd values for each trial, channel, and frequency bin
    psd_speech = zeros(Num_trials, Num_channels, num_freqs);    % psd values for the speech signal for each trial, channel, and frequency bin
    psd_eeg = zeros(Num_trials, Num_channels, num_freqs);       % psd values for the EEG signal for each trial, channel, and frequency bin
    num_speech_samples = length(speech_rawdata.Amplitudes);     % total length of the speech recording, in samples; used to validate that each trial's speech epoch falls entirely within the recording.
%%    for eeg_trial_idx = 2:eeg_struct.Num_trials  
    for eeg_trial_idx = 1:eeg_struct.Num_trials  
        % Speech onset latency is stored in the eeg data structure. 
        speech_onset_latency = eeg_struct.OnsetLatency(eeg_trial_idx);          % onset latency in seconds
        if speech_onset_latency == 0  % this is for the cases when the first onset latency is 0, 
             speech_onset_idx = 1;  % set to 1 to start at the beginning of the speech recording
        else
            speech_onset_idx = round(speech_onset_latency * speech_rawdata.Fs);     % Convert seconds to samples
        end
        speech_offset_idx = speech_onset_idx + speech_epoch_duration_samples-1; % offset is always the same for this dataset. 
        % ------------------------------------------------------------------------------------------------------------
        % Validate that the speech epoch for this trial falls entirely within the bounds of the speech recording.
        % We intentionally do NOT clamp these indices to the valid range: silently clamping would mask a
        % misaligned or out-of-range onset latency by producing a truncated (or even zero/one-sample) speech
        % epoch, which would silently corrupt the PSD/CSD/MSC estimates for that trial (e.g., near-zero PSD from
        % a degenerate epoch leads to 0/0 = NaN in the MSC calculation, which then poisons CSD_chanmeans and
        % MSC_chanmeans for that channel at every frequency, with no indication of why).
        % Instead, we throw a descriptive error so the problem is caught immediately, at the trial that caused it.
        % ------------------------------------------------------------------------------------------------------------
        if speech_onset_idx < 1
            error('Apply2Dataset_CrossSpectralDensity:OnsetBeforeRecordingStart', ...
                ['Trial %d: computed speech onset index (%d samples, %.3f sec) is before the start of the speech recording. ' ...
                 'Check eeg_struct.OnsetLatency for this trial (and any StartTimeOffset/RandomizeOnsets adjustments applied to it).'], ...
                eeg_trial_idx, speech_onset_idx, speech_onset_latency);
        end
        if speech_offset_idx > num_speech_samples
            error('Apply2Dataset_CrossSpectralDensity:OnsetTooCloseToRecordingEnd', ...
                ['Trial %d: speech epoch requires samples %d:%d, but the speech recording is only %d samples ' ...
                 '(%.3f sec) long. Onset latency = %.3f sec, required epoch duration = %.3f sec (%d samples). ' ...
                 'This trial''s onset is too close to (or beyond) the end of the speech recording for the requested epoch duration.'], ...
                eeg_trial_idx, speech_onset_idx, speech_offset_idx, num_speech_samples, ...
                num_speech_samples / speech_rawdata.Fs, speech_onset_latency, ...
                speech_epoch_duration_samples / speech_rawdata.Fs, speech_epoch_duration_samples);
        end

        fprintf("speech onset index: %d, speech offset index: %d, speech epoch duration in samples: %d\n", speech_onset_idx, speech_offset_idx, speech_epoch_duration_samples);
        speech_epoch = speech_rawdata.Amplitudes(speech_onset_idx:speech_offset_idx);  % extract the speech epoch.
        % Prepare the speech for coherence analysis
        % Extract amplitude envelope, bandpass filter, downsample, save phase and magnitude information. 
        fprintf('run preprocessing for speech epoch: trial %d, speech onset latency = %.2f seconds, speech epoch duration = %.2f seconds, number of samples in speech epoch = %d\n', eeg_trial_idx, speech_onset_latency, speech_epoch_duration_seconds, length(speech_epoch));
        [Speech_Struct] = preprocess_speech_epoch(speech_epoch, speech_rawdata.Fs, ...
                                                  highpass_cutoff, lowpass_cutoff, Num_samples_eeg);
        % Loop through EEG channels
        for ch_idx = 1:Num_channels
            eeg_epoch = squeeze(eeg_struct.Data(ch_idx, :, eeg_trial_idx));
            %%%%% figure; plot(eeg_epoch); title(sprintf('Original EEG epoch for trial %d, channel %d', eeg_trial_idx, ch_idx)); pause(0.2);
            % If the eeg_epoch is a row vector, convert it to a column vector for the filtering and CSD calculation functions, which expect column vectors as input.
            % We could do the same thing by setting eeg_epoch = eeg_epoch(:), which will convert it to a column vector regardless of whether it is originally a row or column vector.
            if size(eeg_epoch, 2) ~= 1
                eeg_epoch = eeg_epoch';
            end
            % Band-pass filter EEG epoch
            eeg_epoch_bpf = band_pass_filt(eeg_epoch, eeg_struct.Fs, highpass_cutoff, lowpass_cutoff);
            % Save the cross spectral density values (one value for each analysis frequency) for this trial and channel
            % in a [trials X channels X frequencies] matrix.
            [csd_vals(eeg_trial_idx, ch_idx, :), psd_speech(eeg_trial_idx, ch_idx, :), psd_eeg(eeg_trial_idx, ch_idx, :)] = CrossSpectralDensity(Speech_Struct.envelope, eeg_epoch_bpf, nfft, eeg_struct.Fs);     
        end
    end

    % Save trial-level CSD values in the main results structure.
    CoherenceResults.CSD = csd_vals;
    CoherenceResults.PSD_speech = psd_speech;
    CoherenceResults.PSD_eeg = psd_eeg;

    % Average CSD values across trials at each channel (and each frequency bin) 
    % and save in the main results structure as a [num_channels X num_frequencies] matrix.
    % 
    % Each single-trial CSD value is a complex value, with a magnitude and phase component.  
    % When we average multiple trials together, we'll get an average vector
    % whose magnitude is determined by both 
    % -- the magnitudes of the single trial vectors and
    % -- the similarity (consistency) of their phase values (how much do they point in the same direction).
    % (by phase values, we mean phase differences between the speech and EEG signals at that frequency.) 
    % CSD_chanmeans will be [num_channels X num_freqs]
    CoherenceResults.CSD_chanmeans = squeeze(mean(csd_vals, 1));
    % Average the power spectral density, for the speech and EEG data, across trials, for each channel and frequency.  
    % We will normalize the CSD values by the PSD of the speech and EEG signals at each frequency.
    CoherenceResults.PSD_speech_chanmeans = squeeze(mean(psd_speech, 1));  
    CoherenceResults.PSD_eeg_chanmeans = squeeze(mean(psd_eeg, 1));
    % ------------------------------------------------------------------------------------ 
    % Magnitude Squared Coherence (MSC) provides a measure of coherence 
    % by taking the squared magnitude of the CSD values and normalizing 
    % by the product of the PSD values for the speech and EEG signals at each frequency.
    % MSC = |CSD|^2 / (PSD_speech * PSD_eeg)
    % We'll calculate MSC from CSD and PSD values averaged across trials. 
    % MSC_chanmeans will be [num_channels X num_freqs]
    % ------------------------------------------------------------------------------------ 
    CoherenceResults.MSC_chanmeans  = abs(CoherenceResults.CSD_chanmeans).^2 ./ (CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);                
    % Calculate the frequency vector corresponding to the CSD values and save in the main results structure.
    % This will be the same for all trials and channels.   
    all_freqs = (0:nfft-1) * (eeg_struct.Fs / nfft);   % frequency vector corresponding to all DFT bins; ranges from 0 to fs - fs/nfft
    CoherenceResults.Freqs = all_freqs(1:num_freqs);  % frequencies corresponding to the one-sided spectrum (num_freqs bins, parity-aware -- see above); ranges from 0 to fs/2 for even nfft, or up to just under fs/2 for odd nfft
end
