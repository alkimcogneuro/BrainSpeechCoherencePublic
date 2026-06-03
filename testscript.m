% manage_analysis_cross_spectral_density() is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% it expects a cell array of strings as input, where each string is a filepath to an EEG data file.
% the cell array eegfiles here contains just one filename, which is the EEG data file for subject SE0001 and condition P1.
eegfiles{1} = '~/Downloads/SE0001_P1_targets_expose.mat'

speech_files_path = '~/Downloads/';  % tells us where to find the speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.

%% Our analyses wil focus on frequencies between 2 and 35 Hz, which are the frequencies that are most relevant for speech processing and entrainment.
%% we'll band pass filter the speech envelope between 2 and 35 Hz before computing the cross spectral density, to focus on the frequencies that are most relevant for speech processing and entrainment.
highpass_cutoff = 2.6;  % in Hz. 
lowpass_cutoff = 35;  % in Hz
start_time_offset = 0.35;
epoch_dur = 2;
% This is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);

% if we want to run with optional epoch parameters, we can specify those as follows:
run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  % this will run the CSD analysis using an analysis epoch that starts 0.25 seconds after the original EEG onset and lasts for 2 seconds.
% you don't need to run this load command, unless  you want to inspect the results of the cross spectral density analysis for one subject and one condition,
%  which are saved in a .mat file by the manage_analysis_cross_spectral_density() function.
% This load command loads the results of the cross spectral density analysis for one subject and one condition, which are saved in a .mat file by the manage_analysis_cross_spectral_density() function.
% The loaded variable "Results_CSD_Analysis" is a structure containing the results of the cross spectral density analysis, 
% including the mean CSD values for each channel and frequency, the CSD values for each trial, channel, and frequency, and the frequencies that were analyzed.
% The structure should have the following fields:
% - CSD_chanmeans: a matrix of size [num_channels x num_freqs], where each element is the mean CSD value for that channel and frequency, averaged across all trials.
% - CSD: a 3D matrix of size [num_trials  x num_channels x num_freqs], where each element is the CSD value for that trial, channel, and frequency.
% - Freqs: a vector of length num_freqs, containing the frequencies that were analyzed in the cross spectral density analysis.
% - Chanlocs: a structure containing the channel locations for the EEG data, which can be used for plotting topographic maps of the CSD values.
% - MSC_chanmeans: a matrix of size [num_channels x num_freqs], where each element is the mean MSC value for that channel and frequency, averaged across all trials.
load('~/Downloads/SE0001_P1_targets_expose_cross_spectral_density_results_2_35_Hz.mat');