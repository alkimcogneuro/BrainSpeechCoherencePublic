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
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN   % in seconds
        options.EpochDuration (1,1) double {mustBeReal} = NaN     % in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
    end

    fprintf(' -- Running Cross Spectral Density analysis for subject %s\n', eeg_struct.Subj_id);
    fprintf(' -- Filter parameters: HighPass = %d Hz, LowPass = %d Hz\n', highpass_cutoff, lowpass_cutoff);
    %% fprintf(' --Optional parameters: StartTimeOffset = %.2f seconds, EpochDuration = %.2f seconds\n', options.StartTimeOffset, options.EpochDuration);
    fprintf(' -- EEG sampling rate: %d Hz, Speech sampling rate: %d Hz\n', eeg_struct.Fs, speech_rawdata.Fs);
    fprintf(' -- Number of EEG trials: %d, Number of EEG channels: %d\n', eeg_struct.Num_trials, eeg_struct.Num_channels);
    fprintf(' -- Length of EEG epochs (raw data): %d samples (%.2f seconds)\n', size(eeg_struct.Data, 2), size(eeg_struct.Data, 2) / eeg_struct.Fs);    
    if ~isnan(options.StartTimeOffset)
        fprintf(' --StartTimeOffset: %.2f seconds, which corresponds to %d samples at the EEG sampling rate\n', options.StartTimeOffset, round(options.StartTimeOffset * eeg_struct.Fs));
    else
        fprintf(' --StartTimeOffset: NaN (using full length of EEG epochs)\n');
    end
    if ~isnan(options.EpochDuration)
        fprintf(' --EpochDuration: %.2f seconds, which corresponds to %d samples at the EEG sampling rate\n', options.EpochDuration, round(options.EpochDuration * eeg_struct.Fs));
    else
        fprintf(' --EpochDuration: NaN (using full length of EEG epochs)\n');
    end
    
    Num_channels = size(eeg_struct.Data, 1);  % Check the size of the data matrix; should be [num_channels x num_samples x num_trials]
    Num_samples = size(eeg_struct.Data, 2);
    Num_trials = size(eeg_struct.Data, 3);

    % Z-score normalization of EEG data across all samples and trials, separately for each channel.
    % Reshape to (numchannels x (numsamples*numtrials)) in order to copmute mean and std across all samples and trials for each channel, 
    % then reshape back to original dimensions after normalization.
    eegdata_2d = reshape(eeg_struct.Data, Num_channels, []);
    mu    = mean(eegdata_2d, 2);        % 64x1 mean per channel
    sigma = std(eegdata_2d, 0, 2);      % 64x1 std per channel
    
    % Subtract mean and divide by std, then reshape back to original dimensions
    eeg_data_norm = (eegdata_2d - mu) ./ sigma;
    eeg_struct.Data = reshape(eeg_data_norm, Num_channels, Num_samples, Num_trials);  % reshape back to original dimensions
    % Z-score normalization of speech data.
    % We normalize the speech across the entire speech signal
    % because the speech is not segmented into trials in the same way as the EEG data, and we want to preserve the relative amplitude differences across the entire speech signal,
    speech_amplitudes = speech_rawdata.Amplitudes;
    speech_amplitudes_norm = (speech_amplitudes - mean(speech_amplitudes)) / std(speech_amplitudes);
    speech_rawdata.Amplitudes = speech_amplitudes_norm;  % update the speech raw data structure with the normalized amplitudes.

    % round and convert the start time offset from seconds to samples, 
    % so that we can apply it to the speech onset latencies in the EEG data structure.    
    start_offset_samples = round(options.StartTimeOffset * eeg_struct.Fs); % 
    epoch_duration_samples = round(options.EpochDuration * eeg_struct.Fs); 
    sample_range_restricted = start_offset_samples:(start_offset_samples + epoch_duration_samples - 1);  % the range of samples to use for the analysis epoch, relative to the original EEG onset.  This will be applied to the speech onset latencies to determine the new speech epoch start and end times.

    % If the user provided a StartTimeOffset and EpochDuration, then we will 
    % adjust the EEG data structure to only include the samples in the specified range relative to the original EEG onset.  
    % This will effectively shift the analysis epoch for both the EEG and speech data by the specified offset and duration.  
    % If the user did not provide these optional arguments, 
    % then we will use the full length of the EEG epochs as the analysis epoch duration, 
    % and we will not apply any offset to the speech onset latencies.
    if (~isnan(options.StartTimeOffset) & ~isnan(options.EpochDuration))
        % user provided StartTime — do something
        eeg_struct.Data = eeg_struct.Data(:, sample_range_restricted, :)
        % shift the speech onset latency by the specified offset, so that the new speech epoch will be aligned with the new EEG epoch. 
        % this will work for all trials, because the speech epoch will be defined relative to the new EEG onset for each trial.
        % eeg_struct.OnsetLatency is a vector of length Num_trials, with the original speech onset latency for each trial.
        % we will add the StartTimeOffset (in seconds) to each of these latencies, to get the new speech onset latency for each trial, which will be aligned with the new EEG epoch. 
        eeg_struct.OnsetLatency = eeg_struct.OnsetLatency + options.StartTimeOffset;      
    end
    
    num_samples_epoch = size(eeg_struct.Data, 2);  % the number of samples in the EEG epochs after any optional adjustments to the epoch duration.  This will be used to determine the duration of the speech epochs and to set the nfft parameter for the CSD calculation.
    % nfft is the number of points to use in the FFT calculation for the CSD estimation.
    % we want to set nfft equal to the number of samples in the EEG epochs, so that we are computing the FFT on the full length of the analysis epochs.
    % note that if nfft is a power of 2 for computational efficiency; 
    % this will slightly reduce the frequency resolution but is generally a good idea for FFT calculations.
    % for example, if the EEG epochs are 1500 samples long, then nfft will be set to 1024, which is the largest power of 2 less than or equal to 1500. 
    % This means that the FFT will be computed on the first 1024 samples of the EEG and speech epochs, and the remaining samples will be ignored for the CSD calculation. 
    % nfft = 2^floor(log2(num_samples_epoch));
    nfft = num_samples_epoch;
    % find the largest even number of samples that is less than or equal to num_samples_epoch, to ensure that nfft is even, 
    % which is important for interpreting the frequency bins of the FFT output.
    %if mod(nfft, 2) ~= 0
    %    nfft = nfft - 1;  % make nfft even if it is odd
    %end
    
    fprintf(' --Number of samples in EEG epochs after any optional adjustments: %d samples (%.2f seconds)\n', num_samples_epoch, num_samples_epoch / eeg_struct.Fs);
    fprintf(' --Number of FFT points (nfft) for CSD calculation: %d\n', nfft);  
    CoherenceResults = struct();  % Initialize
    CoherenceResults.analysis_type = 'cross_spectral_density';     
    CoherenceResults.highpass_cutoff = highpass_cutoff;    
    CoherenceResults.lowpass_cutoff = lowpass_cutoff;      
    CoherenceResults.eeg_Fs = eeg_struct.Fs;                               
    CoherenceResults.speech_Fs = speech_rawdata.Fs;                       
    % CoherenceResults. = eeg_struct.eeg_file_name; 
    CoherenceResults.label = eeg_struct.Subj_id;
    CoherenceResults.Chanlocs = eeg_struct.Chanlocs;  

    % Initialize results matrix: NumTrials x NumChannels
    cross_spectral_density_vals = cell(eeg_struct.Num_trials, eeg_struct.Num_channels);
    
    % Duration of speech epoch in seconds, 
    % calculated from the number of samples in the EEG epochs.
    % if the optional arguments for StartTimeOffset and EpochDuration are provided, 
    % then the EEG epoch duration  will have been adjusted in accordance with those parameters 
    % and the speech epoch duration will match that new EEG epoch duration.
    speech_epoch_duration = size(eeg_struct.Data, 2) / eeg_struct.Fs;
    for eeg_trial_idx = 1:eeg_struct.Num_trials
        % Original EEG onset
        % If the user provided a StartTimeOffset, then the EEG onset will have been shifted by that offset, 
        % and the speech epoch will be aligned with the new EEG onset.
        speech_onset_latency = eeg_struct.OnsetLatency(eeg_trial_idx);  
        speech_offset_latency = speech_onset_latency + speech_epoch_duration;

        % Convert seconds to sample indices
        speech_onset_idx = uint64(speech_onset_latency * speech_rawdata.Fs);
        speech_offset_idx = uint64(speech_offset_latency * speech_rawdata.Fs);

        % Ensure indices are within bounds
        speech_onset_idx = max(1, min(speech_onset_idx, length(speech_rawdata.Amplitudes)));
        speech_offset_idx = max(1, min(speech_offset_idx, length(speech_rawdata.Amplitudes)));

        % Extract speech epoch
        % Ren has already cut the speech into epochs that match the EEG epochs, 
        % so this should just be indexing into the correct segment of the speech data.

        speech_epoch = speech_rawdata.Amplitudes(speech_onset_idx:speech_offset_idx);

        % prepare the speech for coherence anaysis
        % Extract amplitude envelope, bandpass filter, downsample, save phase and magnitude information. 

        % i was prevoiusly assuming that eeg_struct would have a filed called Num_samples, but it doesn't,
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
            %% figure; plot(eeg_epoch_bpf); title(sprintf('Trial %d, Channel %d', eeg_trial_idx, ch_idx)); pause(0.2); 
            %% Save the cross spectral density values (one value for each analysis frequency) for this trial and channel
            %% in a trials X channels X frequencies matrix.
            [csd_vals(eeg_trial_idx, ch_idx, :), psd_speech(eeg_trial_idx, ch_idx, :), psd_eeg(eeg_trial_idx, ch_idx, :), mspc(eeg_trial_idx, ch_idx, :)] = CrossSpectralDensity(Speech_Struct.envelope, eeg_epoch_bpf, nfft);     
        end
    end

    % Save trial-level CSD values in the main results structure.
    CoherenceResults.CSD = csd_vals;
    CoherenceResults.PSD_speech = psd_speech;
    CoherenceResults.PSD_eeg = psd_eeg;
    CoherenceResults.MSPC = mspc;

    % AVERAGE ACROSS TRIALS at each channel (and each frequency bin) 
    % and save in the main results structure as a channels X frequencies matrix.
    % 
    % the chanmeans matrices will be num_channels x num_freqs, 
    % where each element is the mean CSD value for that channel and frequency, averaged across all trials.
    % Each single trial CSD value is a complex value, with a magnitude and phase.

    % Note that when we average across trials, consistency of signal-signal phasedifferences across trials
    % will lead to larger average CSD values.
    % If phase is variable, then the vectors will point in different directions across trials, 
    % and the average will be smaller.

    % When we average multiple trials together, we'll get an average vector
    % whose magnitude is determined by the magnitudes of the single trial vectors
    % and the similarity of their phases (how much do they point in the same direction).
    %  
    % the mean of two complex numbers is the mean of their real parts plus i times the mean of their imaginary parts, 
    % so we can just take the mean across trials for each channel and frequency, and the result will be a complex number that represents the average CSD value for that channel and frequency across all trials.
    % this is a vector average, so if the CSD values for a given channel and frequency are consistent across trials (i.e., they have similar phase), 
    % then the average will have a larger magnitude, whereas if the CSD values are variable across trials (i.e., they have different phases), then the average will have a smaller magnitude.
    CoherenceResults.CSD_chanmeans = squeeze(mean(csd_vals, 1));
    CoherenceResults.MSPC_chanmeans = squeeze(mean(mspc, 1));

    % Average the power spectral density, for the speech and EEG data, across trials, for each channel and frequency.  
    % This will allow us to calculate coherence values from the CSD values, 
    % by normalizing the CSD values by the power spectral density of the speech and EEG signals at each frequency.
    CoherenceResults.PSD_speech_chanmeans = squeeze(mean(psd_speech, 1));  
    CoherenceResults.PSD_eeg_chanmeans = squeeze(mean(psd_eeg, 1));
    % We'll calculate the Magnitude Squared Coherence (MSC) values from the CSD and PSD values, averaged across trials, 
    % for each channel and frequency.
    % Calculate coherence from the CSD and PSD values, averaged across trials, for each channel and frequency, 
    % and save in the main results structure as a channels X frequencies matrix.
    % at each channel and frequency: MSC = |CSD|^2 / (PSD_speech * PSD_eeg)

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
    CoherenceResults.MSC_chanmeans  = abs(CoherenceResults.CSD_chanmeans).^2 ./ (CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);    
    % we could calculate "coherence" values by taking the magnitude of the CSD values and normalizing by the PSD values, without squaring the magnitudes,
    % but this would not be the standard definition of coherence, and the resulting values would not be unitless or bounded between 0 and 1, so it would be harder to interpret them as reflecting the strength of the relationship between the signals at each frequency.
    %%    CoherenceResults.MSC_chanmeans  = abs(CoherenceResults.CSD_chanmeans) ./ sqrt(CoherenceResults.PSD_speech_chanmeans .* CoherenceResults.PSD_eeg_chanmeans);
        
    % calculate the frequency vector corresponding to the CSD values, and save in the main results structure.
    % this will be the same for all trials and channels, since we are using the same nfft and sampling rate for all of them.    
    all_freqs = (0:nfft-1) * (eeg_struct.Fs / nfft);  % frequency vector corresponding to DFT bins; ranges from 0 to fs - fs/nfft
    CoherenceResults.Freqs = all_freqs(1:nfft/2 + 1);  % frequencies corresponding to the one-sided spectrum; ranges from 0 to fs/2 
end

% we're calculating CSD_freqs twice, once inside the cSD function and once here....
% CSD_freqs = (0:(nfft/2)) * (1/signal_len);  % frequency vector for the one-sided spectrum, in cycles per sample. 
% the calculation in the CSD function is more accurate, because it takes into account the actual sampling rate of the EEG data, whereas this calculation assumes a sampling rate of 1 Hz.

