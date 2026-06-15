function [] = plot_allchans_MSC_byfreqs(msc_results)    
    % arguments
    % -- msc_results is a structure containing the results of the CSD analysis for one subject and one condition.
    % 
    % The field "MSC_chanmeans" is a [num_channels x num_freqs] matrix, 
    % in which each element is the mean MSC value for one channel x frequency, averaged across trials.
   
    % Calculate the mean MSC at each frequency, averaged across all channels.
    % note: We've already calculated the channel means for MSC at each frequency (when we first calculated MSC).
    msc_avg_channel_means = mean(msc_results.MSC_chanmeans, 1);  % Average across channels (dimension 1)    
    % For each channel, average MSC values within a narrow band of interest.
    msc_lowerbd_freq = 3;   % lower bound of frequency range of interest
    msc_upperbd_freq = 7;   % upper bound of frequency range of interest

    freq_idx = (msc_results.Freqs >= msc_lowerbd_freq) & (msc_results.Freqs <= msc_upperbd_freq);  % logical index for frequencies within your band of interest
    % create a [num_channels X num_freqs] matrix of the absolute values of the MSC_chanmeans, 
    % where num_freqs is the number of frequencies in the narrow band of interest (e.g., 3-7 Hz), 
    % and then average across frequencies within the band of interest (dimension 2), resulting in a vector of length num_channels.
    msc_avg_withinband_by_channel = mean(msc_results.MSC_chanmeans(:, freq_idx), 2);  % average across frequencies (dimension 2), resulting in a vector of length num_channels.

    % We need to convert the values to double precision for plotting, 
    % I think this is because the plotting code assumes that the input values can be very small, 
    % and they will be rounded to zero if we keep them in single precision.
    msc_avg_withinband_by_channel = double(msc_avg_withinband_by_channel);  % convert to double precision for plotting
    fprintf('msc_avg %.2f - %.2f Hz \n', msc_lowerbd_freq, msc_upperbd_freq);
    disp(msc_avg_withinband_by_channel); % print the average MSC values for each channel, averaged across frequencies within the band of interest, to the command window.

    channel_locations = msc_results.Chanlocs;  % grab the channel locations from the CSD results structure, for use in plotting topomap.
    msc_avg_withinband_by_channel = msc_avg_withinband_by_channel(:);   % make sure that we have a column vector, because plot_topomap expects this. 
    % plot_topomap(msc_avg_withinband_by_channel', msc_results.Chanlocs)  % plot the topomap using the coherence values and channel locations from the EEG data file.
    call_topoplot_eeglab(msc_avg_withinband_by_channel', msc_results.Chanlocs, 'CSD', 'CSD')  % plot the topomap using the coherence values and channel locations from the EEG data file.
    % ------------------------------------------------------------------------------ %
    % grab the figure handle for the topoplot that was just created, 
    % and set the title of the figure to indicate the subject ID and condition for this data.
    % save the figure to a file.  
    % ------------------------------------------------------------------------------ %
    fig_handle = gcf;  % get the current figure handle
    % set the title of the figure to indicate the subject ID and condition for this data.
    title(sprintf('MSC %.2f - %.2f Hz, Subject %s, Condition %s, ExposureType %s', msc_lowerbd_freq, msc_upperbd_freq, msc_results.Subj_id, msc_results.Block, msc_results.exposure_type));
    saveas(fig_handle, sprintf('~/Downloads/MSC_topomap_%s_%s_%s_%.2f-%.2fHz.png', msc_results.Subj_id, msc_results.Block, msc_results.exposure_type, msc_lowerbd_freq, msc_upperbd_freq));  % save the figure to a file with a name that includes the subject ID, condition, and frequency range.
    

    for freq_idx = 1:length(msc_results.Chanlocs)
        fprintf('channel %s, msc(2-5Hz) = %d\n',  msc_results.Chanlocs(freq_idx).labels,  msc_avg_withinband_by_channel(freq_idx) );
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
