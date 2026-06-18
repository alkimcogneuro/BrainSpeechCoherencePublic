
function [] = run_CSD_dataset_list(eeg_data_filenames, speech_files_path, highpass_cutoff, lowpass_cutoff, options)
    % This script runs the cross spectral density analysis for all data files in a list of filenames, using the specified filter parameters.
    % For each data file, it runs Apply2Dataset_CrossSpectralDensity() to analyze the entrainment between the EEG data and the corresponding speech data,
    % 
    % input arguments:
    arguments
        eeg_data_filenames (1,:) cell {mustBeText}   % cell array of strings; each string is a filepath to an EEG data file.
        speech_files_path (1,1) string {mustBeText}  % string, path to the folder containing the speech .wav files.  
                                                     % The names of the individual .wav files should match the audio_file field in the EEG data structures.
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
    datetime_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    results_foldername = sprintf('CSD_Analysis_Results_%s', datetime_str);  % folder name for results files
    % grab the path to the first EEG data file in the list, and use that to create a path for the output folder, which will be in the same directory as the EEG data files.
    output_dir = fullfile(fileparts(eeg_data_filenames{1}), results_foldername);  % create a path for the output folder, which is in the same directory as the EEG data files.
    if ~exist(output_dir, 'dir')  % Check if the output folder already exists; it should ont, because we just made up the new folder name using the current date and time, but we check just in case.
        mkdir(output_dir);  % if the output folder does not exist, create it
    end
    %------------------------------------------------------------------------------------------------------%    
    % Iterate through all the EEG data files.
    % Run Apply2Dataset_CrossSpectralDensity() for each one, using the corresponding speech data. 
    % For each EEG data file, we load a data structure called EEG_struct
    % and then we load the corresponding speech data from a .wav file
    % and run Apply2Dataset_CrossSpectralDensity()
    % collect the results into an array of structures called csd_results_structures.
    %------------------------------------------------------------------------------------------------------%    
    for idx = 1:length(eeg_data_filenames)
        fprintf ('Processing eeg data file %d:  %s\n', idx, eeg_data_filenames{idx});
        % fprintf ('output file %d:  %s\n', idx, csd_outputfiles(idx)) 
        eeg_data_file = eeg_data_filenames{idx};    % grab a single eeg data file name
        [eeg_filepath, eeg_filestem, eeg_fileext] = fileparts(eeg_data_file);   % segment the file path into its components
        load(eeg_data_file);                        % Load EEG data structure from one file 
        disp('Size of eeg data structure is: ');  disp(size(EEG_struct.Data));
        % The EEG_struct structure should contain a field called "audio_file" 
        % which names the corresponding speech .wav file.  
        % We use that to load the correct speech data for this subject and condition.
        fprintf('Current EEG Data corresponds to AudioFile: %s\n', EEG_struct.audio_file);
        % Read speech from a file and package it just a bit.
        speech_full_filepath = fullfile(speech_files_path, EEG_struct.audio_file);  % Create full path to speech file
        [speech_amplitudes, fs_speech] = audioread (speech_full_filepath);          % Reaad speech raw data from a .wav file
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
        csd_output_file = sprintf("%s/CSD_results_%s.mat", output_dir, eeg_filestem); 
        % If the user specified optional values for the epoching parameters, we pass those to the function, and if they don't, we call the function without those arguments, 
        % and it will default to the full length of the EEG epochs as the analysis epoch duration.    

        if isnan(options.StartTimeOffset) || isnan(options.EpochDuration)
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, RandomizeOnsets=options.RandomizeOnsets);
        else    
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, ...
             StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration, RandomizeOnsets=options.RandomizeOnsets); 
        end
        % Write each individual subject's result structure, contents of the variable Results_CSD_Analysis, to a file.
        save(csd_output_file, 'Results_CSD_Analysis');
        % We write a sidecar file to accompany each results file, 
        % which records information about the analysis parameters (e.g., filter parameters, epoching parameters, etc.) that we want to keep track of.
        % The sidecar file is in text format for easy viewing and sharing, and it has the same name as the results file but with "_metadata.txt" appended to the end of the filename.
        metadata_output_file_txt = sprintf("%s/CSD_results_Metadata_%s.txt", output_dir, eeg_filestem);  % Create a filename for the metadata text file, which has the same name as the results file but with "_metadata.txt" appended to the end of the filename.
        fid = fopen(metadata_output_file_txt, 'w');  % Open the metadata text file for writing
        fprintf(fid, 'Metadata for Cross Spectral Density results file\nResults filename: %s\n', csd_output_file);  % 
        fprintf(fid, 'EEG data file: %s\n', eeg_data_file);  % Write the name of the EEG data file that this results file corresponds to
        fprintf(fid, 'High-pass filter cutoff frequency: %d Hz\n', highpass_cutoff);  % Write the high-pass filter cutoff frequency used for this analysis
        fprintf(fid, 'Low-pass filter cutoff frequency: %d Hz\n', lowpass_cutoff);    % Write the low-pass filter cutoff frequency used for this analysis
        fprintf(fid, 'Optional Epoch start time offset: %d seconds\n', options.StartTimeOffset);  % Write the epoch start time offset used for this analysis
        fprintf(fid, 'Optional Epoch duration: %d seconds\n', options.EpochDuration);      % Write the epoch duration used for this analysis
        
        fprintf(fid, 'Randomized onsets status: %d\n', options.RandomizeOnsets);  % Write the status of the randomize onsets flag used for this analysis
        fprintf(fid, 'Analysis date and time: %s\n', datetime_str);  % Write the date and time when the analysis was run
        fprintf(fid, 'Analysis script: run_CSD_dataset_list.m\n');  % Write the name of the analysis script used for this analysis
        fprintf(fid, 'Analysis function: Apply2Dataset_CrossSpectralDensity.m\n');  % Write the name of the analysis function used for this analysis
    
        

        fclose(fid);  % Close the metadata text file
    end
end
