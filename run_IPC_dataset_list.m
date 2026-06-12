
function [] = run_IPC_dataset_list(eeg_data_filenames, speech_files_path, highpass_cutoff, lowpass_cutoff, options)
    % input arguments:
    %   eeg_data_filenames: cell array of strings; each string is a filepath to an EEG data file.
    %                       I'm assuming that each EEG data file contains data from one subject and one condition, 
    %                       and that the filename includes the subject ID and condition (e.g., "SE0001_P1_eeg_data.mat").
    %   speech_files_path: string, path to the folder containing the speech .wav files.  
    %                      The names of the individual .wav files should match the audio_file field in the EEG data structures.
    %   highpass_cutoff: scalar, high-pass filter cutoff frequency in Hz (e.g., 2 Hz)
    %   lowpass_cutoff: scalar, low-pass filter cutoff frequency in Hz (e.g., 35 Hz)
    % how to call with optional arguments;
    % 
    arguments
        eeg_data_filenames (1,:) cell {mustBeText}  % cell array of strings; each string is a filepath to an EEG data file.
        speech_files_path (1,1) string {mustBeText}  % string, path to the folder containing the speech .wav files.  
                                     % The names of the individual .wav files should match the audio_file field in the EEG data structures.
        highpass_cutoff (1,1) double {mustBeReal} = 2  % highpass_cutoff defaults to 2 Hz 
        lowpass_cutoff  (1,1) double {mustBeReal} = 35  % lowpass_cutoff defaults to 35 Hz
        options.StartTimeOffset (1,1) double {mustBeReal} = NaN % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.  
                                                                % If NaN, use the original EEG onset as the start of the analysis epoch.
        options.EpochDuration (1,1) double {mustBeReal} = NaN   % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
    end
    
    % We create two equal-length lists of eeg data files and speech datafiles
    % Iterate through the two lists, 
    %    analyzing the entrainment between one eeg data set and one speech dataset.     
    % arguments
    % 
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
    
    %------------------------------------------------------------------------------------------------------%    
    % Iterate through all the EEG data files.
    % And run Apply2Dataset_InstaPhaseCoherence() for each one, using the corresponding speech data. 
    % For each EEG data file, we load a data structure called EEG_struct
    % and then we load the corresponding speech data from a .wav file
    % collect the results into an array of structures called ???.
    %------------------------------------------------------------------------------------------------------%    
    

    % each file records the name of the EEG file that it corresponds to, so we can keep track of which results correspond to which EEG data files.
    % in addition, we can pair the results file with a sidecar file that contains the name of the EEG data file that it corresponds to, 
    % and any other relevant information about the analysis parameters (e.g., filter parameters, epoching parameters, etc.) that we want to keep track of.
    % and these results will be stored inside the results object itself.
    %*****

    % Create a folder for the results of this analysis, 
    % which will hold the results files for analysis of all files in the list eeg_data_filenames
    datetime_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    results_foldername = sprintf('IPC_Analysis_Results_%s', datetime_str);
    % grab the path to the first EEG data file in the list, and use that to create a path for the output folder, which will be in the same directory as the EEG data files.
    output_dir = fullfile(fileparts(eeg_data_filenames{1}), results_foldername);  % create a path for the output folder, which is in the same directory as the EEG data files.
    if ~exist(output_dir, 'dir')  % Check if the output folder already exists; it should ont, because we just made up the new folder name using the current date and time, but we check just in case.
        mkdir(output_dir);  % if the output folder does not exist, create it
    end
    % LOOP THROUGH THE DATA FILES AND ANALYZE EACH.
    for idx = 1:length(eeg_data_filenames)
        fprintf ('Processing eeg data file %d:  %s\n', idx, eeg_data_filenames{idx});
        eeg_data_file = eeg_data_filenames{idx};    % Grab a single eeg data file name
        [eeg_filepath, eeg_filestem, eeg_fileext] = fileparts(eeg_data_file);   % segment the file path into its components
        load(eeg_data_file);                        % Load EEG data structure from one file 
        disp('Size of eeg data structure is: ');  disp(size(EEG_struct.Data));
        % The EEG_struct structure should contain a field called "audio_file" 
        % which names the corresponding speech .wav file.  
        % We use that to load the correct speech data for this subject and condition.
        fprintf('Current EEG Data corresponds to AudioFile: %s\n', EEG_struct.audio_file);
        % Read speech from a file and package it just a bit.
        speech_full_filepath = fullfile(speech_files_path, EEG_struct.audio_file);  % Create full path to speech file
        [speech_amplitudes, fs_speech] = audioread(speech_full_filepath);          % Read speech raw data from a .wav file
        Speech_RawData = struct;                            % Initialize structure
        Speech_RawData.Amplitudes = speech_amplitudes;      % Save the speech amplitudes
        Speech_RawData.Fs = fs_speech;                      % Save the speech sampling rate (44.1 khz for .wav files)
        % Create a filename for saving analysis results to the same folder as the EEG data.
        % The name will be unique to the EEG data filename (e.g., include the subject ID and condition) 
        
        % Each call to Apply2Datset_CrossSpectralDensity returns a structure with the results for one EEG data set and one speech data set.
        % we collect those structures in an array called ipc_results_structures.
        
        % The call to Apply2Dataset_CrossSpectralDensity() will be 
        % contingent on the optional arguments for epoching.  If the user specifies those arguments, we pass them to the function, and if they don't, we call the function without those arguments, and it will use the default values (which are to use the full length of the EEG epochs as the analysis epoch duration, and to use the original EEG onset as the start of the analysis epoch).
        % 
        % The output file will include the name of the folder we created for the results of this analysis, 
        % and the file name will include the the EEG filename (which includes the subject ID and condition)
        % 
        % The filename will NOT include the filter parameters and epoching parameters.  We'll try to keep the filenames simpler. 
        % We will include the meta data in the results structure itself, so that we can keep track of which results correspond to which analysis parameters without having to parse the filename.
        % We will also write a sidecar file that contains the name of the EEG data file that this results file corresponds to, and any other relevant information about the analysis parameters (e.g., filter parameters, epoching parameters, etc.) that we want to keep track of.
        ipc_output_file = sprintf("%s/IPC_results_%s.mat", output_dir, eeg_filestem);
        if isnan(options.StartTimeOffset) || isnan(options.EpochDuration)
            Results_IPC_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff);
        else    
            Results_IPC_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration); 
        end

        % Write each individual subject's result structure, contents of the variable Results_IPC_Analysis, to a file.
        save(ipc_output_file, 'Results_IPC_Analysis');
        % We write a sidecar file (.txt) that contains the name of the EEG data file that this results file corresponds to
        % and other relevant information about the analysis parameters (e.g., filter parameters, epoching parameters, etc.) that we want to keep track of.
        metadata_output_file_txt = sprintf("%s/IPC_results_Metadata_%s.txt", output_dir, eeg_filestem);  % Create a filename for the metadata text file, which has the same name as the results file but with "_metadata.txt" appended to the end of the filename.
        fid = fopen(metadata_output_file_txt, 'w');  % Open the metadata text file for writing
        fprintf(fid, 'Metadata for Instantaneous Phase Coherence analysis results file\nResults filename: %s\n', ipc_output_file);  % 
        fprintf(fid, 'EEG data file: %s\n', eeg_data_file);  % Write the name of the EEG data file that this results file corresponds to
        fprintf(fid, 'High-pass filter cutoff frequency: %d Hz\n', highpass_cutoff);  % Write the high-pass filter cutoff frequency used for this analysis
        fprintf(fid, 'Low-pass filter cutoff frequency: %d Hz\n', lowpass_cutoff);    % Write the low-pass filter cutoff frequency used for this analysis
        fprintf(fid, 'Epoch start time offset: %d seconds\n', options.StartTimeOffset);  % Write the epoch start time offset used for this analysis
        fprintf(fid, 'Epoch duration: %d seconds\n', options.EpochDuration);      % Write the epoch duration used for this analysis
        fclose(fid);  % Close the metadata text file

        
    end
end
