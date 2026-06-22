% manage_analysis_cross_spectral_density() is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% it expects a cell array of strings as input, where each string is a filepath to an EEG data file.
% the cell array eegfiles here contains just one filename, which is the EEG data file for subject '~/DATA/PseudoLanguageFallon/SE0001 and condition P1.
%% eegfiles = {'~/DATA/PseudoLanguageFallon/SE0004_P2_targets_expose.mat', '~/DATA/PseudoLanguageFallon/SE0004_P2_targets_train.mat'};
eegfiles = {
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0001_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0002_E_fullrec.mat', 
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0004_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0005_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0006_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0009_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0010_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0012_E_fullrec.mat'}; 

% skipping subjct 3 because they've got only 63 channels
% skipping Subject 8 because it might have the wrong audio file pointer(eeg epochs included latencies that were beyond the length of the speech recordingt.)
% '~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0008_E_fullrec.mat', 
% '~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0009_E_fullrec.mat', 
% '~/DATA/PseudoLanguageFallon_all_fulldata_files/SE00010_E_fullrec.mat'};
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0011_E_fullrec.mat',

%% eegfiles = {'~/DATA/PseudoLanguageFallon/fake_eeg_data_perfect.mat'};

% cell array of strings; each string is a filepath to an EEG data file. The filenames should include the subject ID and condition (e.g., "SE0001_P1_eeg_data.mat").

speech_files_path = '~/DATA/PseudoLanguageFallon';  % tells us where to find the speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
%% Velveteen_pseudo.wav
%% Our analyses wil focus on frequencies between 2 and 35 Hz, which are the frequencies that are most relevant for speech processing and entrainment.
%% we'll band pass filter the speech envelope between 2 and 35 Hz before computing the cross spectral density, to focus on the frequencies that are most relevant for speech processing and entrainment.
highpass_cutoff = 2.6;  % in Hz. 
lowpass_cutoff = 35;    % in Hz
%% start_time_offset = 0.2;
%% epoch_dur = 2;
% This is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);

% if we want to run with optional epoch parameters, we can specify those as follows:
%%run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  % this will run the CSD analysis using an analysis epoch that starts 0.25 seconds after the original EEG onset and lasts for 2 seconds.

[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);  % This will run the CSD analysis using the full epoch.
[csd_rand_onsets_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, RandomizeOnsets=true);  % run control analysis. 

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


%% load('~/~/DATA/PseudoLanguageFallon/~/DATA/PseudoLanguageFallon/SE0001_P1_targets_expose_cross_spectral_density_results_2_35_Hz.mat');

%{
 
'~/DATA/PseudoLanguageFallon/SE0001_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0001_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0001_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0001_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0002_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0002_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0002_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0002_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0003_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0003_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0003_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0003_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0004_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0004_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0004_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0004_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0005_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0005_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0005_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0005_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0006_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0006_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0006_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0006_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0008_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0008_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0008_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0008_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0009_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0009_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0009_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0009_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0010_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0010_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0010_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0010_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0012_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0012_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0012_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0012_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0013_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0013_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0013_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0013_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0015_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0015_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0015_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0015_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0016_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0016_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0016_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0016_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0017_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0017_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0017_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0017_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0018_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0018_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0018_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0018_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0019_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0019_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0019_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0019_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0020_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0020_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0020_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0020_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0021_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0021_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0021_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0021_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0022_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0022_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0022_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0022_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0023_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0023_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0023_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0023_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0024_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0024_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0024_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0024_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0025_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0025_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0025_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0025_P2_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0026_P1_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0026_P1_targets_train.mat',
'~/DATA/PseudoLanguageFallon/SE0026_P2_targets_expose.mat',
'~/DATA/PseudoLanguageFallon/SE0026_P2_targets_train.mat',
 
%}
