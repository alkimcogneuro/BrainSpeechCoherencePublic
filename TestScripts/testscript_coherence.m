% run_CSD_dataset_list will run the cross spectral density analysis for all data files in the list of files eegfiles.

% It expects a cell array of strings as input, where each string is a filepath to an EEG data file.

% eegfiles ia a cell array of filenames that contain EEG_struct objects
%% eegfiles = {'~/DATA/PseudoLanguageFallon/SE0004_P2_targets_expose.mat', '~/DATA/PseudoLanguageFallon/SE0004_P2_targets_train.mat'};
eegfiles = {
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0001_E_fullrec.mat',
'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0002_E_fullrec.mat'} 
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0004_E_fullrec.mat',
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0005_E_fullrec.mat',
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0006_E_fullrec.mat',
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0009_E_fullrec.mat',
%'~/DATA/PseudoLanguageFallon_all_fulldata_files/SE0010_E_fullrec.mat'}; 

% skipping subjct 3 because they've got only 63 channels
% skipping Subject 8 because it might have the wrong audio file pointer.
% Subject 11 is also problematic.  eeg epochs go beyond the speech.

% eegfiles = {'~/DATA/PseudoLanguageFallon/fake_eeg_data_perfect.mat'};

% cell array of strings; each string is a filepath to an EEG data file. The filenames should include the subject ID and condition (e.g., "SE0001_P1_eeg_data.mat").

speech_files_path = '~/DATA/PseudoLanguageFallon';  % location of the speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
%% Velveteen_pseudo.wav
%% Our analyses wil focus on frequencies between 2 and 35 Hz, which are the frequencies that are most relevant for speech processing and entrainment.
%% we'll band pass filter the speech envelope between 2 and 35 Hz before computing the cross spectral density, to focus on the frequencies that are most relevant for speech processing and entrainment.
highpass_cutoff = 2.6;  % in Hz. 
lowpass_cutoff = 35;    % in Hz
% start_time_offset = 0.2;
% epoch_dur = 2;
% This is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);

% if we want to run with optional epoch parameters, we can specify those as follows:
%%run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  % this will run the CSD analysis using an analysis epoch that starts 0.25 seconds after the original EEG onset and lasts for 2 seconds.

[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);  % This will run the CSD analysis using the full epoch.
%% run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, RandomizeOnsets=true);  % run control analysis. 

