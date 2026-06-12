function [] = rewrite_EEG_struct(filenames)
    %% This function will read a list of filenames, 
    % where each filename is a .mat file containing an EEG_struct structure with the original format.
    % The function will rewrite the EEG_struct structure for each file,
    % to combine the Block and exposure_type fields into a new Condition field,
    % and save the new EEG_struct structure to a new .mat file in a subfolder called "EEG_structs_NEWFORMAT".
    % For example, we could create a list of one EEG structure filenmae
    % with  eegfiles{1} = '~/Downloads/SE0001_P1_targets_expose.mat'
    % and then we would run  
    %  rewrite_EEG_struct(eegfiles)
        
    for idx = 1:length(filenames)
        % re-initialize the new EEG_struct structure for each iteration of the loop 
        % to avoid carrying over any fields from the previous iteration.
        
        filename = filenames{idx};  % grab one filename from the list of filenames
        old_eeg_struct = struct();  % initialize an empty structure to hold the original EEG_struct structure that we will load from the .mat file specified by the filename.
        new_EEG_struct = struct();  % initialize an empty structure to hold the new EEG_struct structure that we will create by copying fields from the original EEG_struct structure and combining the Block and exposure_type fields into a new Condition field.
        old_eeg_struct = load(filename);    % load the EEG_struct structure from the .mat file specified by the filename.
        
        % Grab all the fields from the original EEG_struct structure 
        % and copy them into a new EEG_struct structure.
        new_EEG_struct.OnsetLatency = old_eeg_struct.EEG_struct.OnsetLatency;
        new_EEG_struct.TrigID = old_eeg_struct.EEG_struct.TrigID;
        new_EEG_struct.Trial_cond = old_eeg_struct.EEG_struct.Trial_cond;
        new_EEG_struct.Num_trials = old_eeg_struct.EEG_struct.Num_trials;
        new_EEG_struct.Fs = old_eeg_struct.EEG_struct.Fs;
        new_EEG_struct.Data = old_eeg_struct.EEG_struct.Data;
        new_EEG_struct.Chanlocs = old_eeg_struct.EEG_struct.Chanlocs;
        new_EEG_struct.audio_file = old_eeg_struct.EEG_struct.audio_file;
        new_EEG_struct.Subj_id = old_eeg_struct.EEG_struct.Subj_id;
        new_EEG_struct.Num_channels = old_eeg_struct.EEG_struct.Num_channels;
        %% we don't copy the Block and exposure_type fields from the original EEG_struct structure, 
        %% because we will combine those two fields into a new Condition field in the new EEG_struct structure.
        %% new_EEG_struct.Block = old_eeg_struct.EEG_struct.Block;
        %% new_EEG_struct.exposure_type = old_eeg_struct.EEG_struct.exposure_type;
        
        exposure_type = old_eeg_struct.EEG_struct.exposure_type;    % extract the exposure_type field from the original EEG_struct structure, which indicates whether the condition is "expose" or "unexpose".
        Block = old_eeg_struct.EEG_struct.Block;            % extract the Block field from the original EEG_struct structure, which indicates the experimental block (e.g., "P1", "P2", etc.).
        % Possible values for "Block" are 'E', 'P1', or 'P2' 
        if Block == 'E'
            Block = 'English';  % if the Block field is "E", we keep it as "E" in the new Condition field, since "E" already indicates that this is an exposure condition.
        elseif Block == 'P1'
            Block = 'Pseudowords1';  % if the Block field is "P1", we keep it as "P1" in the new Condition field, since "P1" already indicates that this is a training condition. 
        elseif Block == 'P2'
            Block = 'Pseudowords2';  % if the Block field is "P2", we keep it as "P2" in the new Condition field, since "P2" already indicates that this is a training condition. 
        else
            error('Unexpected value for Block field: %s', Block);  % if the Block field has an unexpected value, we throw an error to alert the user.
        end
        
        % possible values for "exposure_type" are 'e' and 't' (expose/train).
        if exposure_type == 'e'
            exposure_type = 'ExposureNoMeaning';  % if the exposure_type field is "e", we change it to "expose" in the new Condition field, to make it more descriptive and easier to understand.
        elseif exposure_type == 't'
            exposure_type = 'ExposureMeaning';  % if the exposure_type field is "t", we change it to "train" in the new Condition field, to make it more descriptive and easier to understand.
        else
            error('Unexpected value for exposure_type field: %s', exposure_type);  % if the exposure_type field has an unexpected value, we throw an error to alert the user.
        end 
        
        new_cond_field = strcat(Block, '_', exposure_type);  % create a new condition field by combining the Block and exposure_type fields from the original EEG_struct structure.
        new_EEG_struct.Condition = new_cond_field;  % add a new field called "Condition" to the new EEG_struct structure, which combines the Block and exposure_type fields from the original EEG_struct structure.
        
        % Save the old EEG_struct structure to a new .mat file, 
        % with a filename that indicates that this is the original format of the EEG_struct structure, 
        % before we rewrote it to combine the Block and exposure_type fields into a new Condition field. 
        [path, old_filestem, ext] = fileparts(filename);  % extract the filestem from the filename (i.e., the part of the filename without the path or extension)
        % print out values of path, old_filestem, and ext to check that they are correct.
        fprintf('Path: %s\n', path);
        fprintf('Old filestem: %s\n', old_filestem);
        fprintf('Extension: %s\n', ext);
        
        % Save the original EEG_struct structure to a new .mat file 
        % in a subfolder called "backup_old_EEG_structs", 
        % in case we need to revert back to the original format of the EEG_struct structure 
        % for any reason.
        new_path = fullfile(fileparts(filename), 'EEG_structs_NEWFORMAT');  % create a path for the folder, which is in the same directory as the original file.  
        if ~exist(new_path, 'dir')  % check if the new folder already exists
            mkdir(new_path);  % if the new folder does not exist, create it
        end
        % create a full path to the new file by combining the new folder path, the original filestem, and the original extension.
        new_filename = fullfile(new_path, strcat(old_filestem, ext))  % create a filename for new file. Use the same filestem as the original.
        EEG_struct = new_EEG_struct;  % rename structure for saving. When we save the structure, we want its name to be "EEG_struct"
        save(new_filename, 'EEG_struct');
        
    end
end


%{
    We want to change the format of the EEG_struct structure 
    to combine the Block and exposure_type fields into a new Condition field.
    We'll save the new EEG_struct structure to a new .mat file, 

    the original format of the EEG_struct structure is as follows:
 
    EEG_struct = 

  struct with fields:

     OnsetLatency: [39×1 double]
           TrigID: [39×1 double]
       Trial_cond: [39×1 string]
       Num_trials: 39
               Fs: 1000
             Data: [64×3000×39 single]
         Chanlocs: [1×64 struct]
            Block: 'P1'
       audio_file: 'Velveteen_pseudo.wav'
          Subj_id: 'SE0001'
     Num_channels: 64
    exposure_type: "e" 
%}
