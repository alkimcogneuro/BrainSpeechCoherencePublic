

function [] = run_CSD_dataset_list(eeg_data_filenames, speech_files_path, highpass_cutoff, lowpass_cutoff, options)
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
    %
    %   struct with fields:
    %             Data: [63×1500×1688 single]
    %         Chanlocs: [1×63 struct]
    %               Fs: 250
    %       Num_trials: 1688
    %     Num_channels: 63
    %      Num_samples: 1500
    %          Subj_id: 'SE0001'
    %        Condition: 'P1'
    %           TrigID: [1688×1 double]
    %        TrigLabel: [1688×1 string]
    %             Cond: [1688×1 string]
    %     OnsetLatency: [1688×1 double]

    %------------------------------------------------------------------------------------------------------%    
    % Iterate through all the EEG data files.
    % And run Apply2Dataset_CrossSpectralDensity() for each one, using the corresponding speech data. 
    % For each EEG data file, we load a data structure called EEG_struct
    % and then we load the corresponding speech data from a .wav file
    % and run Apply2Dataset_CrossSpectralDensity()
    % collect the results into an array of structures called cross_spectral_density_results_structures.
    %------------------------------------------------------------------------------------------------------%    
    for idx = 1:length(eeg_data_filenames)
        fprintf ('Processing eeg data file %d:  %s\n', idx, eeg_data_filenames{idx});
        % fprintf ('output file %d:  %s\n', idx, cross_spectral_density_outputfiles(idx)) 
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
        % Create a filename for saving analaysis results to the same folder as the EEG data.
        % The name will be unique to the EEG data filename (e.g., include the subject ID and condition) 
        % and should also include the filter parameters used for the analysis, so that we can keep track of which results correspond to which analysis parameters.
        %%%% should also include the epoching parameters, if we add those as optional arguments to Apply2Dataset_CrossSpectralDensity()
        
        cross_spectral_density_output_file = sprintf("%s/%s_cross_spectral_density_results_%d_%d_Hz.mat", eeg_filepath, eeg_filestem, highpass_cutoff, lowpass_cutoff);
        % Each call to Apply2Datset_CrossSpectralDensity returns a structure with the results for one EEG data set and one speech data set.
        % we collect those structures in an array called cross_spectral_density_results_structures.
        
        % The call to Apply2Dataset_CrossSpectralDensity() will be 
        % contingent on the optional arguments for epoching.  If the user specifies those arguments, we pass them to the function, and if they don't, we call the function without those arguments, and it will use the default values (which are to use the full length of the EEG epochs as the analysis epoch duration, and to use the original EEG onset as the start of the analysis epoch).
        if isnan(options.StartTimeOffset) || isnan(options.EpochDuration)
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff);
        else    
            Results_CSD_Analysis = Apply2Dataset_CrossSpectralDensity(EEG_struct, Speech_RawData, highpass_cutoff, lowpass_cutoff, StartTimeOffset=options.StartTimeOffset, EpochDuration=options.EpochDuration); 
        end
        % Write each individual subject's result structure, contents of the variable Results_CSD_Analysis, to a file.
        save(cross_spectral_density_output_file, 'Results_CSD_Analysis');
    end
end
