% script:  runscript_MSC_CodeSwitch.
%
% This script will run the cross spectral density (CSD) analysis for all data files in the list of files eegfiles.
% 
% eegfiles is a cell array of filenames that contain EEG_struct objects
% 
% This script will call run_CSD_dataset_list() with the specified parameters, 
% which will run the CSD analysis for all data files in the specified list eegfiles.

% run_CSD_dataset_list will save its results to a folder with a name that includes the analysis_ID_label, 
% and the current date and time. 
% For example, if analysis_ID_label is "SL_matched_speech", 
% then the results will be saved to a folder named something like 
% "CSD_Analysis_Results_SL_matched_speech_2026_0620_1726_33", 
% where the numbers at the end indicate the date and time when the analysis was run. 
% The results for each subject will be saved to a separate file within that folder, 
% with a name that includes the subject ID and condition .
% (e.g., "CSD_results_CS_struct04.mat")
% 
% We'll set several parameters that will be used by all calls to run_CSD_dataset_list, 
% including the highpass_cutoff and lowpass_cutoff for the band pass filter, 
% and the start_time_offset and epoch_dur for the analysis epochs.

speech_files_path = '~/DATA/CodeSwitch_Data/CriticalAudioFiles';  % location of all speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
highpass_cutoff = 2.5;  % lower bound of band pass filter, in Hz
lowpass_cutoff = 35;    % upper limit of band pass filter, in Hz
start_time_offset = 1.0;    % in seconds. This is the time offset for the start of the analysis epoch, relative to the original EEG epoch onset. 
epoch_dur = 3;          % in seconds. This is the duration of the analysis epoch. The end of the analysis epoch will be determined by adding this duration to the start time offset. So, for example, if start_time_offset is 1.0 and epoch_dur is 3, then the analysis epoch will begin at +1 seconds after the original EEG onset and end at +4 seconds after the original EEG onset.


% eeg files for the SL (same language) condition. 
eegfiles = {
 '~/DATA/CodeSwitch_Data/SLepochs/SL_struct01.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct02.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct04.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct05.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct06.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct07.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct08.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct09.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct10.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct11.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct12.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct13.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct14.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct15.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct16.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct17.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct18.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct19.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct20.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct21.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct22.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct23.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct24.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct25.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct26.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct27.mat', 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct28.mat',
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct29.mat' 
'~/DATA/CodeSwitch_Data/SLepochs/SL_struct30.mat'
} 
analysis_ID_label1 = 'SL_matched_speech';  % label for the main analysis condition, which will be used in the output filenames and figure titles.
analysis_ID_label2 = 'SL_mismatched_speech';  % label for the control analysis condition, which will be used in the output filenames and figure titles. The control analysis will use randomized onsets to break the temporal relationship between the speech and EEG data, while keeping all other parameters the same as the main analysis.
[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label1,  highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  
[csd_control_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label2, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur, RandomizeOnsets=true);  % run control analysis.

% eeg files for the CS (same language) condition. 
eegfiles = {  '~/DATA/CodeSwitch_Data/CSepochs/CS_struct01.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct02.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct04.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct05.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct06.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct07.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct08.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct09.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct10.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct11.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct12.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct13.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct14.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct15.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct16.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct17.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct18.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct19.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct20.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct21.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct22.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct23.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct24.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct25.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct26.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct27.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct28.mat',
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct29.mat', 
    '~/DATA/CodeSwitch_Data/CSepochs/CS_struct30.mat',
} 

analysis_ID_label1 = 'CS_matched_speech';  % label for the main analysis condition, which will be used in the output filenames and figure titles.
analysis_ID_label2 = 'CS_mismatched_speech';  % label for the control analysis condition, which will be used in the output filenames and figure titles. The control analysis will use randomized onsets to break the temporal relationship between the speech and EEG data, while keeping all other parameters the same as the main analysis.
[csd_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label1,  highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur);  
[csd_control_results_files_list] = run_CSD_dataset_list(eegfiles, speech_files_path, analysis_ID_label2, highpass_cutoff, lowpass_cutoff, StartTimeOffset=start_time_offset, EpochDuration=epoch_dur, RandomizeOnsets=true);  % run control analysis.

