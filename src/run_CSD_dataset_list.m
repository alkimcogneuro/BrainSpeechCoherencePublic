

function [result_file_list] = run_CSD_dataset_list(eeg_data_filenames, speech_files_path, analysis_ID_label, highpass_cutoff, lowpass_cutoff, options)
    % This script runs the cross spectral density analysis for all data files in a list of filenames, using the specified filter parameters.
    % For each data file, it runs Apply2Dataset_CrossSpectralDensity() to analyze the entrainment between the EEG data and the corresponding speech data,
    % 
    % input arguments:
    arguments
        eeg_data_filenames (1,:) cell {mustBeText}   % cell array of strings; each string is a filepath to an EEG data file.
        speech_files_path (1,1) string {mustBeText}  % string, path to the folder containing the speech .wav files.  
        % The names of the individual .wav files should match the audio_file field in the EEG data structures.
        analysis_ID_label (1,1) string {mustBeText}  % string label for this analysis, which will be used in the output filenames and figure titles.  This is optional, and defaults to "main_analysis" if not specified.
        highpass_cutoff (1,1) double {mustBeReal} = 2   % highpass_cutoff, defaults to 2 Hz 
        lowpass_cutoff  (1,1) double {mustBeReal} = 35  % lowpass_cutoff, defaults to 35 Hz
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.  
        % If NaN, use the original EEG onset as the start of the analysis epoch.
        options.EpochDuration (1,1) double {mustBeReal} = NaN   % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
        options.RandomizeOnsets (1,1) logical = false  % optional flag for whether to randomize the onset times of the speech epochs relative to the EEG data.  
        % If true, the speech onset times will be randomly changed.  
        % This can be used as a control analysis to test whether any observed entrainment effects are specific to the original alignment of the speech and EEG data, 
        % or if they could be explained by non-specific temporal correlations.
    end
    fprintf('Running run_CSD_dataset_list() for %d EEG data files, with highpass cutoff %.2f Hz and lowpass cutoff %.2f Hz\n', length(eeg_data_filenames), highpass_cutoff, lowpass_cutoff);
    %--------------------------------------------------------------------------------------%
    % EEG data structure will have a format APPROXIMATELY LIKE THIS:
    % EEG data structure will have this format:
    % EEG_struct = 
    %   struct with fields:
    %      OnsetLatency: [39×1 double]
    %            TrigID: [39×1 double]
    %        Trial_cond: [39×1 string]
    %        Num_trials: 39
    %                Fs: 1000
    %              Data: [64×3000×39 single]
    %          Chanlocs: [1×64 struct]
    %             Block: 'P1'               % should change this to "Condition" in the future, to be more general and less tied to our specific experiment.
    %        audio_file: 'Velveteen_pseudo.wav'
    %           Subj_id: 'SE0001'           
    %      Num_channels: 64
    %     exposure_type: "e"
    %--------------------------------------------------------------------------------------%
    % Create a folder for the results of this analysis, 
    % which will hold the results files for analysis of all files in the list eeg_data_filenames
    datetime_str = char(datetime('now', 'Format', 'yyyy_MMdd_HHmm_ss'));
    results_foldername = sprintf('CSD_Analysis_Results_%s_%s', analysis_ID_label, datetime_str);  % folder name for results files
    % grab the path to the first EEG data file in the list, and use that to create a path for the output folder, which will be in the same directory as the EEG data files.
    output_dir = fullfile(fileparts(eeg_data_filenames{1}), results_foldername);  % create a path for the output folder, which is in the same directory as the EEG data files.
    if ~exist(output_dir, 'dir')  % Check if the output folder already exists; it should ont, because we just made up the new folder name using the current date and time, but we check just in case.
        mkdir(output_dir);  % if the output folder does not exist, create it
    end
    fprintf('   Output directory for CSD results: %s\n', output_dir);  % print the output directory to the command window

    %------------------------------------------------------------------------------------------------------%    
    % Iterate through all the EEG data files.
    % Run Apply2Dataset_CrossSpectralDensity() for each one, using the corresponding speech data. 
    % For each EEG data file, we load a data structure called EEG_struct
    % and then we load the corresponding speech data from a .wav file
    % and run Apply2Dataset_CrossSpectralDensity()
    % collect the results into an array of structures called csd_results_structures.
    %------------------------------------------------------------------------------------------------------%    
    result_file_list = [];
    All_Results_CSD_Analysis = struct;  % Initialize a structure to hold the results for all subjects and conditions, which will be saved to a single .mat file at the end of the function.
    for idx = 1:length(eeg_data_filenames)
        fprintf ('  run_CSD_dataset_list() \n\n   Processing eeg data file %d:  %s\n', idx, eeg_data_filenames{idx});
        eeg_data_file = eeg_data_filenames{idx};    % grab a single eeg data file name
        [eeg_filepath, eeg_filestem, eeg_fileext] = fileparts(eeg_data_file);   % segment the file path into its components
        load(eeg_data_file);                        % Load EEG data structure from one file 
        % The EEG_struct structure should contain a field called "audio_file" 
        % which names the corresponding speech .wav file.  
        % We use that to load the correct speech data for this subject and condition.
        fprintf('Current EEG Data corresponds to AudioFile: %s\n', EEG_struct.audio_file);
        % Read speech from a file and package it just a bit.
        speech_full_filepath = fullfile(speech_files_path, EEG_struct.audio_file);  % Create full path to speech file
        fprintf('Audio file full path: %s\n', speech_full_filepath);
        [speech_amplitudes, fs_speech] = audioread(speech_full_filepath);          % Reaad speech raw data from a .wav file
        Speech_RawData = struct;                            % Initialize structure
        Speech_RawData.Amplitudes = speech_amplitudes;      % Save the speech amplitudes
        Speech_RawData.Fs = fs_speech;                      % Save the speech sampling rate (44.1 khz for .wav files)
        % --------------------------------------------------------------------------------------%
        % Each call to Apply2Datset_CrossSpectralDensity returns a structure, which we'll call Results_CSD_Analysis.
        % Results_CSD_Analysis will contain the results of the Cross Spectral Density analysis for one subject and one condition.
        % The Results_CSD_Analysis data structure will include fields for meta data about the analysis parameters,
        % including filter parameters, epoching parameters, etc., so that we can keep track of the analysis parameters that were used for each results file.
        % --------------------------------------------------------------------------------------%
        % We create a filename, csd_output_file for saving the analysis results.  
        % The filename will include the EEG data structure filename, which records the subject ID and some condition information. 
        % --------------------------------------------------------------------------------------%
        csd_output_file = fullfile(output_dir, sprintf("CSD_results_%s.mat", eeg_filestem));   
        result_file_list = [result_file_list, csd_output_file];  % collect results in this list, which will be returned as the output of this function.
            
        % If the user specified optional values for the epoching parameters, we pass those to the function, and if they don't, we call the function without those arguments, 
        % and it will default to the full length of the EEG epochs as the analysis epoch duration.    

        Results_CSD_Analysis = struct;  % Initialize the results structure for this subject and condition, which will be filled with the output of Apply2Dataset_CrossSpectralDensity() and then saved to a file.
        if isnan(options.StartTimeOffset) || isnan(options.EpochDuration)
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, RandomizeOnsets=options.RandomizeOnsets);
        else    
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, ...
            StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration, RandomizeOnsets=options.RandomizeOnsets); 
        end
        % Write each individual subject's result structure, contents of the variable Results_CSD_Analysis, to a file.
        save(csd_output_file, 'Results_CSD_Analysis');
        % In addition to saving the individual subjects analyses to seprarate files, 
        % we will also collect all the results into a structure array called All_Results_CSD_Analysis, 
        % where each element of the array corresponds to one subject and condition, 
        % and has a field called Results_CSD_Analysis which contains the results for that subject and condition.
        % we'll save the All_Results_CSD_Analysis structure to a single .mat file at the end of the function, 
        % so that we have all the results in one place for easy loading and sharing.
        All_Results_CSD_Analysis(idx).Results_CSD_Analysis = Results_CSD_Analysis;  % Save the results for this subject and condition into the All_Results_CSD_Analysis structure, which will be saved to a file at the end of the function.
    end
    
    all_results_file = fullfile(output_dir, "All_Results_CSD_Analysis.mat");  % Create a filename for the All_Results_CSD_Analysis structure, which will be saved to a .mat file in the output directory.
    save(all_results_file, 'All_Results_CSD_Analysis');    % Save the All_Results_CSD_Analysis structure to a .mat file in the output directory.
            
    % Save meta data for this run.
    % this will include:
    %   -- the list of results files
    %   -- the list of input EEG data files.
    %   -- the filter parameters used for this analysis (e.g., high-pass and low-pass cutoff frequencies)
    %   -- the epoching parameters used for this analysis (e.g., epoch start time offset and epoch duration)
    %   -- the date and time when the analysis was run, which will be useful for keeping track of different runs of the analysis,
    %  -- the name of the analysis script and function used for this analysis, which will be useful for keeping track of the code that was used to generate the results,
    % and any other parameters that we want to keep track of for this analysis.  
    metadata_output_file = fullfile(output_dir, "CSD_results_Metadata.txt");  % Create a filename for the metadata text file for this run
    fid = fopen(metadata_output_file, 'w');  % Open the metadata text file for writing
    fprintf(fid, 'Metadata for Cross Spectral Density (and Magnitude Squared Coherence) Analysis\n\n');  % 
    fprintf(fid, 'High-pass filter cutoff frequency: %.2f Hz\n', highpass_cutoff);  % Write the high-pass filter cutoff frequency used for this analysis
    fprintf(fid, 'Low-pass filter cutoff frequency: %.2f Hz\n', lowpass_cutoff);    % Write the low-pass filter cutoff frequency used for this analysis
    fprintf(fid, 'Optional Epoch start time offset: %d seconds\n', options.StartTimeOffset);  % Write the epoch start time offset used for this analysis
    fprintf(fid, 'Optional Epoch duration: %d seconds\n', options.EpochDuration);      % Write the epoch duration used for this analysis
    fprintf(fid, 'Randomized onsets status: %d\n', options.RandomizeOnsets);  % Write the status of the randomize onsets flag used for this analysis
    fprintf(fid, 'Analysis date and time: %s\n', datetime_str);  % Write the date and time when the analysis was run
    fprintf(fid, 'Analysis script: run_CSD_dataset_list.m\n');  % Write the name of the analysis script used for this analysis
    fprintf(fid, 'Analysis function: Apply2Dataset_CrossSpectralDensity.m\n');  % Write the name of the analysis function used for this analysis    

    % write the list of eeg files that were analyzed to a text file, which will be useful for keeping track of the results files.
    fprintf(fid, '\n--List of EEG data files analyzed:\n');  % Write a header for the list of EEG data files that were analyzed
    for eegidx = 1:length(eeg_data_filenames)
        eeg_data_file = eeg_data_filenames{eegidx};  % grab the name of the EEG data file for this index
        fprintf(fid, '%s\n', eeg_data_file);  % Write the name of the EEG data file that this results file corresponds to
    end

    % write the list of results files that were generated to the text file, which will be useful for keeping track of the results files
    fprintf(fid, '--List of CSD results files generated:\n');  % Write a header for the list of results files that were generated
    for result_idx = 1:length(result_file_list)
       fprintf(fid, '%s\n', result_file_list{result_idx});  % Write the name of the EEG data file that this results file corresponds to
    end
    fprintf(fid, 'All results files saved to: %s\n', all_results_file);  % Print the name of the file that contains all the results for this run of the analysis
    fclose(fid);  % Close the metadata text file

end
