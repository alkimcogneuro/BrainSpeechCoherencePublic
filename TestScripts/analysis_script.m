

% grab all results structures.
results_file_CS = '~/DATA/CodeSwitch_Data/CSepochs/CSD_Analysis_Results_CS_matched_speech_2026_0620_1729_58/All_Results_CSD_Analysis.mat';
results_file_SL = '~/DATA/CodeSwitch_Data/SLepochs/CSD_Analysis_Results_SL_matched_speech_2026_0620_1726_33/All_Results_CSD_Analysis.mat';

results_file_CS_mismatch = '~/DATA/CodeSwitch_Data/CSepochs/CSD_Analysis_Results_CS_mismatched_speech_2026_0620_1731_30/All_Results_CSD_Analysis.mat';
results_file_SL_mismatch = '~/DATA/CodeSwitch_Data/SLepochs/CSD_Analysis_Results_SL_mismatched_speech_2026_0620_1728_11/All_Results_CSD_Analysis.mat';

% Each load produces an "all_results" structure, which has a field All_Results_CSD_Analysis, which is a 1 x num_subjects struct array, where each element is the results for one subject. Each subject's results have a field Results_CSD_Analysis, which contains the fields CSD_chanmeans, PSD_speech_chanmeans, PSD_eeg_chanmeans, Freqs, and Chanlocs, among others.
all_results_CS = load(results_file_CS);
all_results_SL = load(results_file_SL);
all_results_CS_mismatch = load(results_file_CS_mismatch);
all_results_SL_mismatch = load(results_file_SL_mismatch);

band_low = 3;
band_high = 7;
Chanlocs = all_results_CS.All_Results_CSD_Analysis(1).Results_CSD_Analysis.Chanlocs;  % grab channel locations from the first subject's results, for use in plotting topography. We assume that all subjects have the same channels in the same order, so we can just use the Chanlocs from the first subject.
num_subjects = length(all_results_CS.All_Results_CSD_Analysis);

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


for idx = 1:num_subjects
    fprintf('Subject %d: %s\n %s\n', idx, all_results_CS.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id, all_results_SL.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id);
    subject_IDs_CS{idx} = all_results_CS.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id;
    [MSC_band_CS, band_freqs, band_bin_idx] = BandAverageMSC(all_results_CS.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_CS(idx, :) = MSC_band_CS;    % collect the band-averaged MSC values for this subject and condition into the group matrix, where each row is a subject and each column is a channel.
    [MSC_band_SL, band_freqs, band_bin_idx] = BandAverageMSC(all_results_SL.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_SL(idx, :) = MSC_band_SL;    % collect the band-averaged MSC values for this subject and condition into the group matrix, where each row is a subject and each column is a channel.

    [MSC_band_CS_mismatch, band_freqs, band_bin_idx] = BandAverageMSC(all_results_CS_mismatch.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_CS_mismatch(idx, :) = MSC_band_CS_mismatch; 
    [MSC_band_SL_mismatch, band_freqs, band_bin_idx] = BandAverageMSC(all_results_SL_mismatch.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_SL_mismatch(idx, :) = MSC_band_SL_mismatch;
end 


PlotGroupMSC_Topography(MSC_band_group_CS, MSC_band_group_SL, Chanlocs, Condition1_Label="After Code Switch", Condition2_Label="Same Language Control", ... 
    PlotIndividualSubjects=true, SubjectIDs=subject_IDs_CS, IndividualOutputDir='~/Downloads/MSC_topomaps_by_subject_CS/');
PlotGroupMSC_Topography(MSC_band_group_CS, MSC_band_group_CS_mismatch, Chanlocs, Condition1_Label="After Code Switch", Condition2_Label="After Code Switch, Mismatching Speech Control");
PlotGroupMSC_Topography(MSC_band_group_SL, MSC_band_group_SL_mismatch, Chanlocs, Condition1_Label="Same Language", Condition2_Label="Same Language Mismatching Speech Control");


%{
 % grab all results structures.
results_file_match = '~/DATA/CodeSwitch_Data/SLepochs/CSD_Analysis_Results_matched_20260619_110053/All_Results_CSD_Analysis.mat';
results_file_control = '~/DATA/CodeSwitch_Data/SLepochs/CSD_Analysis_Results_mismatched_20260619_110116/All_Results_CSD_Analysis.mat';

% each load produces an "all_results" structure, which has a field All_Results_CSD_Analysis, which is a 1 x num_subjects struct array, where each element is the results for one subject. Each subject's results have a field Results_CSD_Analysis, which contains the fields CSD_chanmeans, PSD_speech_chanmeans, PSD_eeg_chanmeans, Freqs, and Chanlocs, among others.
all_results_match = load(results_file_match);
all_results_control = load(results_file_control);

band_low = 3.5;
band_high = 7;

Chanlocs = all_results_match.All_Results_CSD_Analysis(1).Results_CSD_Analysis.Chanlocs;  % grab channel locations from the first subject's results, for use in plotting topography. We assume that all subjects have the same channels in the same order, so we can just use the Chanlocs from the first subject.
num_subjects = length(all_results_match.All_Results_CSD_Analysis);
MSC_band_group_match = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the main analysis condition.
MSC_band_group_control = zeros(num_subjects, length(Chanlocs));  % preallocate a matrix to hold the band-averaged MSC values for each subject and each channel, for the control analysis condition.

for idx = 1:num_subjects
    fprintf('Subject %d: %s\n %s\n', idx, all_results_match.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id, all_results_control.All_Results_CSD_Analysis(idx).Results_CSD_Analysis.Subj_id);

    [MSC_band_match, band_freqs, band_bin_idx] = BandAverageMSC(all_results_match.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_match(idx, :) = MSC_band_match;

    [MSC_band_control, band_freqs, band_bin_idx] = BandAverageMSC(all_results_control.All_Results_CSD_Analysis(idx).Results_CSD_Analysis, band_low, band_high)
    MSC_band_group_control(idx, :) = MSC_band_control;
end 

PlotGroupMSC_Topography(MSC_band_group_match, MSC_band_group_control, Chanlocs, Condition1_Label="matched", Condition2_Label="mismatched");
 
%}
