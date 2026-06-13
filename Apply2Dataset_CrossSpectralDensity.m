function [CoherenceResults] = Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, options)
    % ===============================================================================================================
    % This function applies Cross Spectral Density analysis to one EEG dataset. 
    %
    % Arguments:
    %   eeg_struct:     structure containing EEG data and metadata
    %   speech_rawdata: structure containing raw speech data and metadata
    %   highpass_cutoff: high-pass filter cutoff frequency (Hz)
    %   lowpass_cutoff:  low-pass filter cutoff frequency (Hz)
    %
    % example call with optional arguments provided
    % Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff, StartTimeOffset=0.2, EpochDuration=1.8)
    % Returns:
    %   CoherenceResults: structure containing results
    % 
    % functionality to add:
    % think about how to flexibly adjust the width of the analysis epochs for speech and eeg
    % we may need to use sub-sections of the epochs that are provided. 
    % ===============================================================================================================
    arguments
        eeg_struct
        speech_rawdata
        highpass_cutoff (1,1) double {mustBeReal} = 2
        lowpass_cutoff  (1,1) double {mustBeReal} = 35
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.  
                                                                % If NaN, use the original EEG onset as the start of the analysis epoch.

        options.EpochDuration (1,1) double {mustBeReal} = NaN   % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
    end

    fprintf(' -- Running Cross Spectral Density analysis for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Filter parameters: HighPass = %d Hz, LowPass = %d Hz\n', highpass_cutoff, lowpass_cutoff);
    %% fprintf(' --Optional parameters: StartTimeOffset = %.2f seconds, EpochDuration = %.2f seconds\n', options.StartTimeOffset, options.EpochDuration);
    fprintf(' -- EEG sampling rate: %d Hz, Speech sampling rate: %d Hz\n', eeg_struct.Fs, speech_rawdata.Fs);
    fprintf(' -- Number of EEG trials: %d, Number of EEG channels: %d\n', eeg_struct.Num_trials, eeg_struct.Num_channels);
    fprintf(' -- Length of EEG epochs (raw data): %d samples (%.2f seconds)\n', size(eeg_struct.Data, 2), size(eeg_struct.Data, 2) / eeg_struct.Fs);    
    
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
    CoherenceResults.analysis_type = 'cross_spectral_density';     
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
    for eeg_trial_idx = 1:eeg_struct.Num_trials
        % Speech onset latency is stored in the eeg data structure. 
        speech_onset_latency = eeg_struct.OnsetLatency(eeg_trial_idx);  
        speech_offset_latency = speech_onset_latency + speech_epoch_duration;
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
        [Speech_Struct] = preprocess_speech_epoch(speech_epoch, speech_rawdata.Fs, ...
                                                  highpass_cutoff, lowpass_cutoff, num_samples);
        % Loop through EEG channels
        for ch_idx = 1:eeg_struct.Num_channels
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
            % in a trials X channels X frequencies matrix.
            [csd_vals(eeg_trial_idx, ch_idx, :), psd_speech(eeg_trial_idx, ch_idx, :), psd_eeg(eeg_trial_idx, ch_idx, :), mspc(eeg_trial_idx, ch_idx, :)] = CrossSpectralDensity(Speech_Struct.envelope, eeg_epoch_bpf, nfft);     
        end
    end

    % Save trial-level CSD values in the main results structure.
    CoherenceResults.CSD = csd_vals;
    CoherenceResults.PSD_speech = psd_speech;
    CoherenceResults.PSD_eeg = psd_eeg;
    CoherenceResults.MSPC = mspc;

    % Average CSD values across trials at each channel (and each frequency bin) 
    % and save in the main results structure as a [num_channels X num_frequencies] matrix.
    % 
    % Each single-trial CSD value is a complex value, with a magnitude and phase component.  
    % When we average multiple trials together, we'll get an average vector
    % whose magnitude is determined by both 
    % -- the magnitudes of the single trial vectors and
    % -- the similarity (consistency) of their phase values (how much do they point in the same direction).
    % (by phase values, we mean phase differences between the speech and EEG signals at that frequency.) 
    CoherenceResults.CSD_chanmeans = squeeze(mean(csd_vals, 1));
    CoherenceResults.MSPC_chanmeans = squeeze(mean(mspc, 1));
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
    % ------------------------------------------------------------------------------------ 
    CoherenceResults.MSC_chanmeans  = abs(CoherenceResults.CSD_chanmeans).^2 ./ (CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);            
    % Calculate the frequency vector corresponding to the CSD values and save in the main results structure.
    % This will be the same for all trials and channels.   
    all_freqs = (0:nfft-1) * (eeg_struct.Fs / nfft);   % frequency vector corresponding to all DFT bins; ranges from 0 to fs - fs/nfft
    CoherenceResults.Freqs = all_freqs(1:nfft/2 + 1);  % frequencies corresponding to the one-sided spectrum; ranges from 0 to fs/2 
end

% we're calculating CSD_freqs twice, once inside the cSD function and once here....



    % The MSC values will reflect the strength and consistency of the relationship between the speech and EEG signals 
    % at each frequency,
    % The MSC values are the CSD values normalized by the power spectral density of the speech and EEG signals at each frequency.
    % The MSC values are influenced by the magnitudes of the CSD values (the numerator of the MSC calculation), 
    %   and also influenced by the product of the PSD values (the denominator of the MSC calculation),
    % The MSC values will be higher when the CSD values are larger, 
    %   but they will also be higher when the product of the PSD values is smaller, 
    %   which can occur when the power of the signals at that frequency is low.
    % The MSC values will be unitless and will range from 0 to 1, 
    % with 0 indicating no relationship between the signals at that frequency, 
    % and 1 indicating a perfect relationship between the signals at that frequency.
    % 
    % Remember: 
    % The CSD values themselves (MSC numerator) can be influenced by the overall power of the signals at each frequency, 
    % as well as the consistency of the phase relationship between the signals across trials, 
    % which means that they can be higher when the power of the signals at that frequency is high, 
    % and/or when the phase relationship between the signals at that frequency is consistent across trials.
    % 
    % The CSD values and PSD values are averaged across trials, 
    % so the MSC values will reflect the consistency of phase-differences across trials, 
    % as well as the strength of the relationship between the speech and EEG signals at each frequency.
    % we do not normalize the CSD values by the PSD values before averaging across trials, 
    % because we want the trial-level CSD values to reflect the consistency of the power and phase relationships across trials, 
    % and then we can calculate the MSC values from the averaged CSD and PSD values, 
    % which will reflect the overall strength and consistency of the relationship between the speech and EEG signals at each frequency.
    % The squared magnitude of the CSD values is used in the numerator because 
    % coherence is defined as the squared magnitude of the cross-spectral density 
    % normalized by the product of the power spectral densities of the two signals.
    % the product of the power spectral densities is used in the denominator 
    % because it normalizes the CSD values by the overall power of the signals at each frequency,
    % so that the coherence values reflect the strength of the relationship between the signals at each frequency,
    % note that we squared magnitudes in the numerator, resulting in power values, 
    % but we did not square the PSD values in the denominator, so the units of the MSC values will be different from the units of the CSD values.
    %----
    % Note about the unitless values:
    % The CSD magnitudes are in units of power (since we squared the magnitudes), 
    % and the PSD values are also in units of power, 
    % so when we square the CSD magnitudes and divide by the product of the PSD values,
    % the resulting MSC values will be unitless, and will range from 0 to 1,
