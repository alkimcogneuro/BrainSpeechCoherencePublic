% it expects a cell array of strings as input, where each string is a filepath to an EEG data file.
% the cell array eegfiles here contains just one filename, which is the EEG data file for subject '~/DATA/PseudoLanguageFallon/SE0001 and condition P1.
eegfiles = {  '~/DATA/PseudoLanguageFallon/SE0001_P1_targets_expose.mat', 
        '~/DATA/PseudoLanguageFallon/SE0002_P1_targets_train.mat'};

speech_files_path = '~/DATA/PseudoLanguageFallon';  % tells us where to find the speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
%% Velveteen_pseudo.wav
%% Our analyses wil focus on frequencies between 2 and 35 Hz, which are the frequencies that are most relevant for speech processing and entrainment.
%% we'll band pass filter the speech envelope between 2 and 35 Hz before computing the cross spectral density, to focus on the frequencies that are most relevant for speech processing and entrainment.
highpass_cutoff = 2.6;  % in Hz. 
lowpass_cutoff = 35;    % in Hz
%%start_time_offset = 0.2;
%%epoch_dur = 2;
% This is the main function that runs the cross spectral density analysis for all data files in the specified folder, using the specified filter parameters.
% run_CSD_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);

% if we want to run with optional epoch parameters, we can specify those as follows:

%run_IPC_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  % this will run the CSD analysis using an analysis epoch that starts 0.25 seconds after the original EEG onset and lasts for 2 seconds.

run_IPC_dataset_list(eegfiles, speech_files_path, highpass_cutoff, lowpass_cutoff);  % this will run the CSD analysis using the full epoch.

%{
   '~/DATA/PseudoLanguageFallon/SE0002_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0002_P2_targets_train.mat',
  '~/DATA/PseudoLanguageFallon/SE0003_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0003_P2_targets_train.mat'
  '~/DATA/PseudoLanguageFallon/SE0004_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0004_P2_targets_train.mat',
  '~/DATA/PseudoLanguageFallon/SE0005_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0005_P2_targets_train.mat',
  '~/DATA/PseudoLanguageFallon/SE0006_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0006_P2_targets_train.mat',


  '~/DATA/PseudoLanguageFallon/SE0008_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0008_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0009_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0009_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0010_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0010_P2_targets_train.mat',


  '~/DATA/PseudoLanguageFallon/SE0012_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0012_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0013_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0013_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0015_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0015_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0016_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0016_P2_targets_train.mat',

  '~/DATA/PseudoLanguageFallon/SE0017_P2_targets_expose.mat',
  '~/DATA/PseudoLanguageFallon/SE0017_P2_targets_train.mat', 
%}