function [] = plot_allchans_CSD_byfreqs(csd_results)    
    % csd_results is a structure containing the results of the cross spectral density analysis for one subject and one condition.
    % the field "CSD_chanmeans" is a matrix of size [num_channels x num_freqs], where each element is the mean CSD value for that channel and frequency, averaged across all trials.
    % the field "CSD" is a 3D matrix of size [num_trials  x num_channels x num_freqs], 
    % where each element is the CSD value for that trial, channel, and frequency
    %

    % csd_results.CSD is num_trials x num_channels x num_freqs, 
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
    csd_avg_mags = squeeze(mean(abs(csd_results.CSD), [1 2])); 
    size(csd_avg_mags)

    % we already calculated the channel means for CSD at each frequency (when we first calculated CSD).
    % now, calculate the mean CSD at each frequency, averaged across all channels.
    csd_avg_channel_means = abs(mean(csd_results.CSD_chanmeans, 1));  % average across channels (dimension 1)    

    % for each channel, average CSD magnitudes between 2 and 5 Hz and average them.
    csd_lowerbd_freq = 3;  % lower bound of frequency range of interest
    csd_upperbd_freq = 8;  % upper bound of frequency range of interest

    idx = (csd_results.Freqs >= csd_lowerbd_freq) & (csd_results.Freqs <= csd_upperbd_freq);  % logical index for frequencies between 2 and 5 Hz

    csd_avg_2_5Hz_by_channel = mean(abs(csd_results.CSD_chanmeans(:, idx)), 2);  % average across frequencies (dimension 2), resulting in a vector of length num_channels.
    % we need to convert to double precision for plotting, 
    % I think this is because the plotting code assumes that the input values can be very small, 
    % and when we plot them, they will be rounded to zero if we keep them in single precision.
    csd_avg_2_5Hz_by_channel = double(csd_avg_2_5Hz_by_channel);  % convert to double precision for plotting
    fprintf('csd_avg 2-5 Hz \n');
    channel_locations = csd_results.Chanlocs; 
    csd_avg_2_5Hz_by_channel = csd_avg_2_5Hz_by_channel(:);   % this conversion might not be necessary, but it ensures that the input to plot_topomap is a column vector, which is what the function expects.
    %% plot_topomap(csd_avg_2_5Hz_by_channel', csd_results.Chanlocs)  % plot the topomap using the random voltages and the channel locations from the EEG data file.
    call_topoplot_eeglab(csd_avg_2_5Hz_by_channel', csd_results.Chanlocs,'CSD', 'CSD')  % plot the topomap using the random voltages and the channel locations from the EEG data file.
    for idx = 1:length(csd_results.Chanlocs)
        fprintf('channel %s, CSD(2-5Hz) = %d\n',  csd_results.Chanlocs(idx).labels,  csd_avg_2_5Hz_by_channel(idx) );
    end

    % plot the MSC values (which are the coherence values) for each channel, averaged across frequencies between 2 and 5 Hz, in a topographic map.
    % these are based on CSD averaged across trials.
    msc_avg_2_5Hz_by_channel = mean(abs(csd_results.MSC_chanmeans(:, idx)), 2);  % average across frequencies (dimension 2), resulting in a vector of length num_channels.
    % we need to convert to double precision for plotting, 
    % I think this is because the plotting code assumes that the input values can be very small, 
    % and when we plot them, they will be rounded to zero if we keep them in single precision.
    msc_avg_2_5Hz_by_channel = double(msc_avg_2_5Hz_by_channel);  % convert to double precision for plotting
    fprintf('csd_avg 2-5 Hz \n');
    channel_locations = csd_results.Chanlocs; 
    msc_avg_2_5Hz_by_channel = msc_avg_2_5Hz_by_channel(:);   % this conversion might not be necessary, but it ensures that the input to plot_topomap is a column vector, which is what the function expects.

    call_topoplot_eeglab(msc_avg_2_5Hz_by_channel', csd_results.Chanlocs,'MSC', 'MSC')  % plot the topomap using the random voltages and the channel locations from the EEG data file.
    for idx = 1:length(csd_results.Chanlocs)
        fprintf('channel %s, MSC(2-5Hz) = %d\n',  csd_results.Chanlocs(idx).labels,  msc_avg_2_5Hz_by_channel(idx) );
    end

    figure; 
    for channel = 1:length(csd_results.Chanlocs) 
        chanlabel = csd_results.Chanlocs(channel).labels;
        %%Results_CSD_Analysis.Chanlocs.labels
        csd_channel_means = csd_results.CSD_chanmeans(channel, :);
        csd_channel_mags = abs(csd_channel_means);
        freqs = csd_results.Freqs;
        % Identify the indices of frequencies between 2 and 20 Hz, 
        % so that we can plot the CSD values in that range.
        idx = (freqs >= 0.5) & (freqs <= 20);
        subplot(ceil(length(csd_results.Chanlocs)/4), 4, channel);
        % grab the magnitude of the CSD values for the specified channel, and plot them against the frequencies:
        plot(freqs(idx), csd_channel_mags(idx), 'LineWidth', 2); 
        xlabel('Frequency (Hz)');
        ylabel('Cross-Spectral Density');
        title(sprintf('Cross-Spectral Density for Channel %d %s', channel, chanlabel));
    end

    % figure; 
    % for channel = 1:length(csd_results.Chanlocs) 
    %     chanlabel = csd_results.Chanlocs(channel).labels;
    %     %%Results_CSD_Analysis.Chanlocs.labels
    %     msc_channel_means = csd_results.MSC_chanmeans(channel, :);
        
    %     freqs = csd_results.Freqs;
    %     % Identify the indices of frequencies between 2 and 20 Hz, 
    %     % so that we can plot the CSD values in that range.
    %     idx = (freqs >= 0.5) & (freqs <= 20);
    %     subplot(ceil(length(csd_results.Chanlocs)/4), 4, channel);
    %     % grab the magnitude of the CSD values for the specified channel, and plot them against the frequencies:
    %     plot(freqs(idx), msc_channel_means(idx), 'LineWidth', 2); 
    %     xlabel('Frequency (Hz)');
    %     ylabel('Magnitude Squared Coherence');
    %     title(sprintf('MSC for Channel %d %s', channel, chanlabel));
    % end

end
