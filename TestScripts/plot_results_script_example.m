

% This script is for plotting the results of the CSD analysis for the Code Switch dataset. 
% It will load the results structures for all subjects and conditions, 
% extract the band-averaged MSC values for each subject and each channel, 
% and then plot the group-level topography of MSC values for each condition, as well as the difference between conditions. 

% grab all results structures.   
% Note:  user must make sure that the results files are from the correct results folder
% for instance, the folder CSD_Analysis_Results_CS_matched_speech_2026_0620_1729_58
% contains the results for the main analysis of the CS condition, with matched speech controls, 
% and was created by a run on 2026-06-20 at 17:29:58. 

results_file_CS = '~/DATA/MyExperiment_Data/CSepochs/CSD_Analysis_Results_CS_matched_speech_2026_0622_2133_27/All_Results_CSD_Analysis.mat';
results_file_CS_mismatch = '~/DATA/MyExperiment_Data/CSepochs/CSD_Analysis_Results_CS_mismatched_speech_2026_0622_2134_42/All_Results_CSD_Analysis.mat';

results_file_SL = '~/DATA/MyExperiment_Data/SLepochs/CSD_Analysis_Results_SL_matched_speech_2026_0622_2130_53/All_Results_CSD_Analysis.mat';
results_file_SL_mismatch = '~/DATA/MyExperiment_Data/SLepochs/CSD_Analysis_Results_SL_mismatched_speech_2026_0622_2132_10/All_Results_CSD_Analysis.mat';

% Each loaded file contains an 'All_Results_CSD_Analysis' structure, which is a 1 x num_subjects struct array, 
% where each element is the results for one subject. 
% Each subject's results have a field Results_CSD_Analysis, 
% which contains the fields CSD_chanmeans, PSD_speech_chanmeans, PSD_eeg_chanmeans, Freqs, and Chanlocs, among others.
S = load(results_file_CS);
all_results_CS = S.All_Results_CSD_Analysis

S = load(results_file_SL);
all_results_SL = S.All_Results_CSD_Analysis;

S = load(results_file_CS_mismatch);
all_results_CS_mismatch = S.All_Results_CSD_Analysis;

S = load(results_file_SL_mismatch);
all_results_SL_mismatch = S.All_Results_CSD_Analysis;

band_low = 3;     % what frequency band do we want to average over?  For example, 3-8 Hz is the theta band, which is relevant for speech processing and entrainment.
band_high = 12;
Chanlocs = all_results_CS(1).Results_CSD_Analysis.Chanlocs;  % grab channel locations from the first subject's results, for use in plotting topography. We assume that all subjects have the same channels in the same order, so we can just use the Chanlocs from the first subject.
num_subjects = length(all_results_CS);

MSC_band_group_CS = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the main analysis condition.
MSC_band_group_SL = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the SL analysis condition.
MSC_band_group_CS_mismatch = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the main analysis condition.
MSC_band_group_SL_mismatch = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the SL analysis condition.


% Calculate MSC within a specific frequency band for each subject and each condition, and store the results in the group matrices. 
% Each row of the group matrices corresponds to one subject, and each column corresponds to one channel.

% create a cell array of strings to hold the subject IDs for each subject and each condition, which we'll use for printing individual subject data.
% we want to extract the information in this field: all_results_CS.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id
subject_IDs_CS = cell(num_subjects, 1);
subject_IDs_SL = cell(num_subjects, 1);

% Loop through each subject and extract the band-averaged MSC values for each condition, and store them in the group matrices.
for idx = 1:num_subjects
    fprintf('Subject %d: %s\n %s\n', idx, all_results_CS(idx).Results_CSD_Analysis.Subj_id, all_results_SL(idx).Results_CSD_Analysis.Subj_id);
   
    subject_IDs_CS{idx} = all_results_CS(idx).Results_CSD_Analysis.Subj_id;

    [MSC_band_CS, band_freqs, band_bin_idx] = BandAverageMSC(all_results_CS(idx).Results_CSD_Analysis, band_low, band_high, AvgMethod='power-weighted')
    MSC_band_group_CS(idx, :) = MSC_band_CS;    % collect the band-averaged MSC values for this subject and condition into the group matrix, where each row is a subject and each column is a channel.
    [MSC_band_SL, band_freqs, band_bin_idx] = BandAverageMSC(all_results_SL(idx).Results_CSD_Analysis, band_low, band_high, AvgMethod='power-weighted')
    MSC_band_group_SL(idx, :) = MSC_band_SL;    % collect the band-averaged MSC values for this subject and condition into the group matrix, where each row is a subject and each column is a channel.

    [MSC_band_CS_mismatch, band_freqs, band_bin_idx] = BandAverageMSC(all_results_CS_mismatch(idx).Results_CSD_Analysis, band_low, band_high,  AvgMethod='power-weighted')
    MSC_band_group_CS_mismatch(idx, :) = MSC_band_CS_mismatch; 
    [MSC_band_SL_mismatch, band_freqs, band_bin_idx] = BandAverageMSC(all_results_SL_mismatch(idx).Results_CSD_Analysis, band_low, band_high,  AvgMethod='power-weighted')
    MSC_band_group_SL_mismatch(idx, :) = MSC_band_SL_mismatch;

    subj_ids{idx} = all_results_CS(idx).Results_CSD_Analysis.Subj_id
end 

% PlotGroupMSC_Topography is a function that takes the group matrices of band-averaged MSC values for two conditions, 
% as well as the channel locations, 
% and plots the topography of MSC values for each condition, as well as the difference between conditions.
% PlotGroupMSC_Topography plots MSC values that are averaged across subjects, for each channel.
PlotGroupMSC_Topography(MSC_band_group_CS, MSC_band_group_SL, Chanlocs, Condition1_Label_Label="After Code Switch", Condition2_Label="Same Language Control");
PlotGroupMSC_Topography(MSC_band_group_CS, MSC_band_group_CS_mismatch, Chanlocs, Condition1_Label="After Code Switch", Condition2_Label="After Code Switch, Mismatching Speech Control");
PlotGroupMSC_Topography(MSC_band_group_SL, MSC_band_group_SL_mismatch, Chanlocs, Condition1_Label="Same Language", Condition2_Label="Same Language Mismatching Speech Control");

% Plot the individual subject topographies for each condition, and save the figures to the specified output directory.
%{
 PlotIndividualSubjectsMSC_Topography(MSC_band_group_CS, MSC_band_group_SL, Chanlocs, Condition1_Label="After Code Switch", Condition2_Label="Same Language Control", SubjectIDs=subj_ids, OutputDir='~/Downloads/MSC_Topography_IndividualSubjects_CS_vs_SL');
PlotIndividualSubjectsMSC_Topography(MSC_band_group_CS, MSC_band_group_CS_mismatch, Chanlocs, Condition1_Label="After Code Switch", Condition2_Label="After Code Switch, Mismatching Speech Control", SubjectIDs=subj_ids, OutputDir='~/Downloads/MSC_Topography_IndividualSubjects_CS_vs_CS_mismatch'); 
PlotIndividualSubjectsMSC_Topography(MSC_band_group_SL, MSC_band_group_SL_mismatch, Chanlocs, Condition1_Label="Same Language", Condition2_Label="Same Language Mismatching Speech Control", SubjectIDs=subj_ids, OutputDir='~/Downloads/MSC_Topography_IndividualSubjects_SL_vs_SL_mismatch'); 
 
%}

