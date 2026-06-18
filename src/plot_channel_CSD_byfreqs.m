function [] = plot_channel_CSD_byfreqs(csd_results, channel)    
    % csd_results is a structure containing the results of the cross spectral density analysis for one subject and one condition.
    % the field "cross_spectral_density_chanmeans" is a matrix of size [num_channels x num_freqs], where each element is the mean CSD value for that channel and frequency, averaged across all trials.
    % the field "cross_spectral_density_vals" is a 3D matrix of size [num_trials  x num_channels x num_freqs], 
    % where each element is the CSD value for that trial, channel, and frequency
    %%  Results_CSD_Analysis.Chanlocs(:).labels  
    % i want to calculate the average magnitude of the CSD values, averaged across all channels and trials, for each frequency.
    csd_avg_mags = squeeze(mean(abs(csd_results.cross_spectral_density_vals), [1 2]));  % take magnitude first, then average across trials and channels
    % the result is a vector of length num_freqs, where each element is the average magnitude of the CSD values at that frequency, averaged across all channels and trials.
    size(csd_avg_mags)

    for chidx = 1:64 
        label = Results_CSD_Analysis.Chanlocs(chidx).labels
        channel_labels{string(label)} = chidx;
    end

    csd_channel_means = csd_results.cross_spectral_density_chanmeans(channel, :);
    csd_channel_mags = abs(csd_channel_means);

    % Calculate the mean CSD at each frequency, averaged across all channels.
    csd_avg_channel_means = abs(mean(csd_results.cross_spectral_density_chanmeans, 1));  % average across channels (dimension 1)    

    freqs = csd_results.CSD_freqs;
    % Identify the indices of frequencies between 2 and 20 Hz, 
    % so that we can plot the CSD values in that range.
    
    idx = (freqs >= 0.5) & (freqs <= 20);
    figure;
    % grab the magnitude of the CSD values for the specified channel, and plot them against the frequencies:
    plot(freqs(idx), csd_channel_mags(idx), 'LineWidth', 2); 
    xlabel('Frequency (Hz)');
    ylabel('Cross-Spectral Density');
    title(sprintf('Cross-Spectral Density for Channel %d', channel));
    hold on; 
    % plot the average CSD across all channels on the same plot, for comparison:
    plot(freqs(idx), csd_avg_mags(idx), 'LineWidth', 2, 'Color', 'r');
end
