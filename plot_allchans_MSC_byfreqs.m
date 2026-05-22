function [] = plot_allchans_MSC_byfreqs(msc_results)    
    % msc_results is a structure containing the results of the cross spectral density analysis for one subject and one condition.
    % the field "CSD_chanmeans" is a matrix of size [num_channels x num_freqs], where each element is the mean msc value for that channel and frequency, averaged across all trials.
    % the field "msc" is a 3D matrix of size [num_trials  x num_channels x num_freqs], 
    % where each element is the msc value for that trial, channel, and frequency
    %

    % msc_results.CSD is num_trials x num_channels x num_freqs, 
    % where each element is the CSD value for that trial, channel, and frequency.
    % Calculate the average magnitude of the CSD values, averaged across all channels and trials, for each frequency.
    % Take magnitude first, then average across trials and channels
    % The result is a vector of length num_freqs, 
    % where each element is the average magnitude of the CSD values at that frequency, 
    % averaged across all channels and trials.
    % Note that  we are taking the magnitude first...we're getting rid of phase information.
    % if CSD is bunch of pages,
    % this will calculate the average magnitude of the CSD values across all channels and trials, 
    % for each frequency, resulting in a vector of length num_freqs. 
   
    % we already calculated the channel means for msc at each frequency (when we first calculated msc).
    % now, calculate the mean msc at each frequency, averaged across all channels.
    msc_avg_channel_means = mean(msc_results.MSC_chanmeans, 1);  % average across channels (dimension 1)    

    % for each channel, average MSC values between 2 and 5 Hz and average them.
    msc_lowerbd_freq = 4;  % lower bound of frequency range of interest
    msc_upperbd_freq = 10;  % upper bound of frequency range of interest

    freq_idx = (msc_results.Freqs >= msc_lowerbd_freq) & (msc_results.Freqs <= msc_upperbd_freq);  % logical index for frequencies between 2 and 5 Hz

    %    msc_avg_2_5Hz_by_channel = mean(abs(msc_results.MSC_chanmeans(:, freq_idx)), 2);  % average across frequencies (dimension 2), resulting in a vector of length num_channels.
    msc_avg_2_5Hz_by_channel = mean(msc_results.MSC_chanmeans(:, freq_idx), 2);  % average across frequencies (dimension 2), resulting in a vector of length num_channels.

    % we need to convert to double precision for plotting, 
    % I think this is because the plotting code assumes that the input values can be very small, 
    % and when we plot them, they will be rounded to zero if we keep them in single precision.

    msc_avg_2_5Hz_by_channel = double(msc_avg_2_5Hz_by_channel);  % convert to double precision for plotting
    fprintf('msc_avg 2-5 Hz \n');
    class(msc_avg_2_5Hz_by_channel)  
    disp(msc_avg_2_5Hz_by_channel);

    channel_locations = msc_results.Chanlocs; 
    msc_avg_2_5Hz_by_channel = msc_avg_2_5Hz_by_channel(:);   % this conversion might not be necessary, but it ensures that the input to plot_topomap is a column vector, which is what the function expects.
    %% plot_topomap(msc_avg_2_5Hz_by_channel', msc_results.Chanlocs)  % plot the topomap using the coherence values and channel locations from the EEG data file.
    call_topoplot_eeglab(msc_avg_2_5Hz_by_channel', msc_results.Chanlocs)  % plot the topomap using the coherence values and channel locations from the EEG data file.

    for freq_idx = 1:length(msc_results.Chanlocs)
        fprintf('channel %s, msc(2-5Hz) = %d\n',  msc_results.Chanlocs(freq_idx).labels,  msc_avg_2_5Hz_by_channel(freq_idx) );
    end


    figure; 
    for channel = 1:length(msc_results.Chanlocs) 
        chanlabel = msc_results.Chanlocs(channel).labels;
        %%Results_CSD_Analysis.Chanlocs.labels
        msc_channel_means = msc_results.MSC_chanmeans(channel, :);
        
        freqs = msc_results.Freqs;
        % Identify the indices of frequencies between 2 and 20 Hz, 
        % so that we can plot the CSD values in that range.
        
        idx = (freqs >= 0.5) & (freqs <= 20);
        subplot(ceil(length(msc_results.Chanlocs)/4), 4, channel);
        % grab the magnitude of the CSD values for the specified channel, and plot them against the frequencies:
        plot(freqs(idx), msc_channel_means(idx), 'LineWidth', 2); 
        xlabel('Frequency (Hz)');
        ylabel('Magnitude Square Coherence');
        title(sprintf('Magnitude Squared Coherence for Channel %d %s', channel, chanlabel));
    end 


end
