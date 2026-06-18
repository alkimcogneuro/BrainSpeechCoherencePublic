
% ---------------------------------------------------------- % 
% This code is for playing Jacqui's speech stimuli, so we can check that the triggers
% are in the right place -- that the sounds are playing at the correct latencies relative to the EEG data.
% We read the CSV file into a table
% The CSV file has headers in the first row, 
% and the columns are "set", "word_target", "latency_1", ... "latency_6"
% where "set" indicates the audio recording
% There are additional columns but they are not relevant for this analysis.
% ---------------------------------------------------------- % 

% CSV file containing the triggers for the speech stimuli, including the latencies for each repetition of each word target.
% edit this filename, to match the location of the CSV file on your computer. The CSV file should have headers in the first row, and the columns should be "set", "word_target", "latency_1", ... "latency_6", where "set" indicates the audio recording, and the latency columns indicate the latencies for each repetition of each word target.

triggers_file = '~/DATA/PseudoLanguageFallon/train_test_triggers.csv';  

speech_stim_trigs = readtable(triggers_file);
% iterate through the rows of the table of word_target triggers, 
% and for each row, we check the "set" column to determine which audio file to play, 
% and then we check the latency columns to determine when to play the audio segments.
for i = 12:height(speech_stim_trigs)
    if strcmp(speech_stim_trigs.set{i}, 'vr_1') || strcmp(speech_stim_trigs.set{i}, 'vr_2')  % check if the "set" column for the current row is "vr_1" or "vr_2", which indicates that this row corresponds to the first or second audio recording. 
        % check if the "set" column for the current row is "vr_1", which indicates that this row corresponds to the first audio recording.
        audiofile = '~/DATA/PseudoLanguageFallon/Velveteen_pseudo.wav'  % if the "set" column is "vr_1", we set the audio file to be the Velveteen_pseudo.wav file, which is the audio recording for the first set of stimuli.
    elseif strcmp(speech_stim_trigs.set{i}, 'oz_1') || strcmp(speech_stim_trigs.set{i}, 'oz_2')  % check if the "set" column for the current row is "vr_3", which indicates that this row corresponds to the second audio recording.
        audiofile = '~/DATA/PseudoLanguageFallon/Oz_pseudo.wav';
    else
        fprintf("Row %d: audiofile code is not recognized, skipping playback.\n", i);  % if the "set" column for the current row has an unexpected value, we print a message to indicate that we are skipping playback for this row.
        continue;  % Skip to the next iteration of the loop, since we don't know which audio file to play for this row.
    end
    fprintf("\n-----\nRow %d: word=  %s  , audiofile=%s (%s)\n", i, speech_stim_trigs.word_target{i}, speech_stim_trigs.set{i}, audiofile); 


    for rep = 1:6
        latency_colname = sprintf('latency_%d', rep);  % construct the column name for the current repetition (e.g., "latency_1", "latency_2", etc.)
        speech_onset_time = speech_stim_trigs.(latency_colname)(i);  % access the latency value for the current row and repetition using dynamic field names
        if ~isnan(speech_onset_time)  % check if the latency value is not NaN (i.e., it is a valid number)
            fprintf("  Repetition %d: latency = %.2f seconds\n", rep, speech_onset_time);
            pause(0.5);  % pause for a short duration (e.g., 0.5 seconds) before playing the audio segment, to avoid the print statement to be read
            play_speech_segment(audiofile, speech_onset_time);  % if the latency value is valid, play the corresponding segment of the speech audio file using the play_speech_segment function.
        else
            fprintf("  Repetition %d: latency value is NaN, skipping playback.\n", rep);
        end
    end
end
 
