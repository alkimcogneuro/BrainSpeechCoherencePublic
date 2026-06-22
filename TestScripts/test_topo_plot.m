% Load one of Jacqui's EEG data files (we assume the file contains a data structure called 'EEG_struct'), 
% We just need to grab the Chanlocs field and save it for later use.

% we can save the channel locations information as a separate .mat file, so that we can load it into the topoplot function in EEGLAB.
% set the file path to point to one of your data files. 
load('/Users/aakim/BrainSpeechEntrainment/EEG_FullRecord_Data/SE0003_E_fullrec.mat');
% save the channel locations as a separate .mat file, so that we can load it into the topoplot function in EEGLAB.
channel_locations = EEG_struct.Chanlocs
save('~/Downloads/chan_locs_struct.mat', 'channel_locations')
load('~/Downloads/chan_locs_struct.mat');   % we can load the channel lcoations back in whenever we need them.

voltages1 = zeros(1,64);        % initialize all voltages to 0.
voltages1(29:31) = 1;  % set C3, C1, Cz have a voltage of 1.
voltages1(61:63) = -1;  % set some of the occipital channels to -1, just for kicks.
plot_topomap(voltages1, channel_locations)  % plot the topomap using the random voltages and the channel locations from the EEG data file.

% let's see what a random montage of voltages looks like: 
voltages_rnd = randn(1,63);  % create a random vector of voltages for each channel, to use as input to the topoplot function.
plot_topomap(voltages_rnd, channel_locations)  % plot the topomap using the random voltages and the channel locations from the EEG data file.
