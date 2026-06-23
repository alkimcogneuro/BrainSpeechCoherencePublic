function [] = plot_speech_powerspectra()
    fprintf('CALCULATING AVERAGE POWER SPECTRUM FOR THE Jacqui DATA\n');
    % Read speech files from disk.
    % Create a matrix of speech envelope vectors.
    speech_files = ["/Users/aakim/BrainSpeechEntrainment/SpeechData/Oz_English.wav", ...
    "/Users/aakim/BrainSpeechEntrainment/SpeechData/Oz_pseudo.wav", ...
    "/Users/aakim/BrainSpeechEntrainment/SpeechData/Velveteen_English.wav", ...
    "/Users/aakim/BrainSpeechEntrainment/SpeechData/Velveteen_pseudo.wav"];
    
    for idx = 1:numel(speech_files) 
        fprintf('reading %dth speech file, %s\n', idx, speech_files(idx));
        [speech_data(idx).amplitudes, speech_data(idx).fs] = audioread(speech_files(idx));     % loads speech raw data from a .wav file
        speech_powerspectra(idx).spectrum = calculate_speech_powerspectrum(speech_data(idx).amplitudes, speech_data(idx).fs);       
    end
    
    for idx = 1:numel(speech_powerspectra)
        %%speech_mean_powerspectrum = mean(speech_powerspectra, 1);  % average across rows (speech files)
        curr_spectrum = speech_powerspectra(idx).spectrum;
        fs = speech_data(idx).fs;
        N = size(curr_spectrum, 1); % number of samples
        freqs = (fs/N)*(0:floor((N-1)/2));    % Reference frequencies in Hz.
        idx_vlf_start = floor(1/(fs/N));
        idx_vlf_end = floor(15/(fs/N));
        idx_speech_freq_start = floor(50/(fs/N));
        idx_speech_freq_end = floor(4000/(fs/N));
        idx_wide_band_start  = floor(50/(fs/N));
        idx_wide_band_end = floor(10000/(fs/N));
        figure; 
        fh = tiledlayout(3,1);
        title(fh,  sprintf('power spectrum for file %s', speech_files(idx)));
        ph = nexttile; stem(freqs(idx_vlf_start:idx_vlf_end), curr_spectrum(idx_vlf_start:idx_vlf_end)); 
        title(ph, 'Average Power spectrum of the speech epochs, very low frequency');    xlabel(ph, 'Frequency (Hz)'); ylabel('Power');
        ph = nexttile; stem(freqs(idx_speech_freq_start:idx_speech_freq_end), curr_spectrum(idx_speech_freq_start:idx_speech_freq_end)); 
        title(ph, 'Average Power spectrum of the speech epochs, speech band');  xlabel(ph, 'Frequency (Hz)'); ylabel('Power');
        ph = nexttile; stem(freqs(idx_wide_band_start:idx_wide_band_end), curr_spectrum(idx_wide_band_start:idx_wide_band_end)); 
        title(ph, 'Average Power spectrum of the speech epochs, speech band');  xlabel(ph, 'Frequency (Hz)'); ylabel('Power');
    end
end
