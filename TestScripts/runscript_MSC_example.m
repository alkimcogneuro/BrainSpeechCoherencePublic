% script:  runscript_MSC_MyStudy.
%
% This script will run the cross spectral density (CSD) analysis for all data files in the list of files eegfiles.
% 
% eegfiles is a cell array of filenames that contain EEG_struct objects
% 
% This script will call run_CSD_dataset_list() with the specified parameters, 
% which will run the CSD analysis for all data files in the specified list eegfiles.

% run_CSD_dataset_list will save its results to a folder with a name that includes the analysis_ID_label, 
% and the current date and time. 
% For example, if analysis_ID_label is "Condition1_matched_speech", 
% then the results will be saved to a folder named something like 
% "CSD_Analysis_Results_Condition1_matched_speech_2026_0620_1726_33", 
% where the numbers at the end indicate the date and time when the analysis was run. 
% The results for each subject will be saved to a separate file within that folder, 
% with a name that includes the subject ID and condition .
% (e.g., "CSD_results_Condition2_struct04.mat")
% 
% We'll set several parameters that will be used by all calls to run_CSD_dataset_list, 
% including the highpass_cutoff and lowpass_cutoff for the band pass filter, 
% and the start_time_offset and epoch_dur for the analysis epochs.

speech_files_path = '~/DATA/MyStudy_Data/CriticalAudioFiles';  % location of all speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
highpass_cutoff = 2.5;  % lower bound of band pass filter, in Hz
lowpass_cutoff = 35;    % upper limit of band pass filter, in Hz
start_time_offset = 1.0;    % in seconds. This is the time offset for the start of the analysis epoch, relative to the original EEG epoch onset. 
epoch_dur = 3;          % in seconds. This is the duration of the analysis epoch. The end of the analysis epoch will be determined by adding this duration to the start time offset. So, for example, if start_time_offset is 1.0 and epoch_dur is 3, then the analysis epoch will begin at +1 seconds after the original EEG onset and end at +4 seconds after the original EEG onset.


% eeg files for condition=Condition1. 
eegfiles = {
 '~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct01.mat', 
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct02.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct04.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct05.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct06.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct07.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct08.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct09.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_struct10.mat',
} 
analysis_ID_label1 = 'Condition1_matched_speech';  % label for the main analysis condition, which will be used in the output filenames and figure titles.
analysis_ID_label2 = 'Condition1_mismatched_speech';  % label for the control analysis condition, which will be used in the output filenames and figure titles. The control analysis will use randomized onsets to break the temporal relationship between the speech and EEG data, while keeping all other parameters the same as the main analysis.
[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label1,  highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  
[csd_control_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label2, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur, RandomizeOnsets=true);  % run control analysis.

% eeg files for condition=Condition2. 
eegfiles = {  '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct01.mat', 
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct02.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct04.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct05.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct06.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct07.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct08.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct09.mat',
    '~/DATA/MyStudy_Data/Condition2epochs/Condition2_struct10.mat',
} 

analysis_ID_label1 = 'Condition2_matched_speech';  % label for the main analysis condition, which will be used in the output filenames and figure titles.
analysis_ID_label2 = 'Condition2_mismatched_speech';  % label for the control analysis condition, which will be used in the output filenames and figure titles. The control analysis will use randomized onsets to break the temporal relationship between the speech and EEG data, while keeping all other parameters the same as the main analysis.
[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label1,  highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  
[csd_control_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label2, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur, RandomizeOnsets=true);  % run control analysis.

