
function [] = play_speech_segment(audiofile, onset_time)
    % arguments:
    % -- audiofile:  a filename for a .wav audiofile (e.g., '~/DATA/PseudoLanguageFallon/Velveteen_pseudo.wav')
    % -- onset_time: a latency in seconds (e.g., 1.5), which indicates when to play the audio segment relative to the start of the audio file.

    [audio_data, Fs] = audioread(audiofile);
    % Convert the onset time from seconds to samples
    onset_sample = round(onset_time * Fs);
    % Define the duration of the segment to play (e.g., 2 seconds)
    segment_duration = 2; % in seconds
    segment_samples = round(segment_duration * Fs);
    % Extract the segment of audio data to play
    segment_data = audio_data(onset_sample:(onset_sample + segment_samples - 1));
    % Play the audio segment
    sound(segment_data, Fs);
    % pause for the duration of the segment to allow it to finish playing before the function ends
    pause(segment_duration+0.5);  % add a short buffer to ensure the audio finishes playing
end

