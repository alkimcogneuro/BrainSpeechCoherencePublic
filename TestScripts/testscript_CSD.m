%% --- Parameters for sample input signals ---

%% --- Generate example signals ---
% x is a sum of two sinusoids
% y is a phase-shifted version of x, plus noise
% This gives us a known ground truth to verify the CSD against

f1 = 3;    % Frequency of component 1 (Hz)
f2 = 15;   % Frequency of component 2 (Hz)
phi = pi/4; % Phase offset between x and y at f1

fs_eeg = 1000;          % sampling frequency (Hz)
T = 3;              % duration of signal (seconds)
num_samples_eeg = T * fs_eeg;        % Total number of samples
t_eeg = (0:num_samples_eeg-1) / fs_eeg;  % Time vector

fs_speech = 44100;
num_samples_speech = T*fs_speech;
t_speech = (0:num_samples_speech-1) / fs_speech;

% y shares both frequency components with x, but the f1 component
% is phase-shifted by pi/4. 

% eeg = sin(2*pi*f1*t_eeg + phi) + 0.5 * sin(2*pi*f2*t_eeg) + 0.5*randn(1, num_samples_eeg);
% speech = sin(2*pi*f1*t_speech) + 0.5 * sin(2*pi*f2*t_speech) + 0.5*randn(1, num_samples_speech);

speech = sin(2*pi*f1*t_speech) + 0.5 * sin(2*pi*f2*t_speech);
eeg = sin(2*pi*f1*t_eeg + phi) + 0.5 * sin(2*pi*f2*t_eeg);

speech_analytic_signal = hilbert(speech);
speech_envelope = abs(speech_analytic_signal);


% i want to plot the spectra for the speech, speech envelope, and EEG signals to verify that they contain the expected frequency components at 3 and 15 Hz, 
% and to check for any unexpected peaks that might be introduced by the resampling process.
% i want to focus on the lower frequencies, below 40 Hz, since that's where the expected frequency components are and where we would expect to see any artifacts from the resampling process.
% please write that code below.  now. 
figure
subplot(3,1,1);
pwelch(speech, [], [], 0:0.1:40, fs_speech);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Original Speech Signal');
subplot(3,1,2);
pwelch(speech_envelope, [], [], 0:0.1:40, fs_speech);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Speech Envelope');
subplot(3,1,3);
pwelch(eeg, [], [], 0:0.1:40, fs_eeg);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Simulated EEG Signal');



% Design linear-phase FIR lowpass filter
fir_order = 2048;    % high order = sharp rolloff
f_cutoff = 30;   % cutoff frequency (Hz)
b = fir1(fir_order, f_cutoff / (fs_speech/2), 'low', kaiser(fir_order+1, 8));
% Apply zero-phase filter
%%speech_envelope = filtfilt(b, 1, speech_envelope);

speech2 = sin(2*pi*f1*t_eeg);
speech2_envelope = abs(hilbert(speech2));

% Normalize the speech envelope
% By dividing the deviations from the mean by the standard deviation, we are scaling the envelope so that it has a standard deviation of 1. 
%% speech_envelope = (speech_envelope - mean(speech_envelope)) / std(speech_envelope);
% i think that normalizing the speech envelope is introducing negative values into the envelope. is that true?
% if (any(~isfinite(speech_envelope)))
%     % this could happen if there's a division by zero.
%     error('ERROR . \nSpeech envelope contains non-finite values\n');
% end

% highpass_cutoff = 4;   % Hz
% lowpass_cutoff = 8;   % Hz
% speech_envelope_bpf = band_pass_filt(speech_envelope, fs_speech, highpass_cutoff, lowpass_cutoff);  % filter

% Downsample from speech sampling rate (e.g,. 44.1 kHz) to the EEG sample rate (e.g., 1000 Hz).
% so we can compare the speech envelope and EEG vectors directly in the CSD analysis.

% unlike our phase coherence analysis, we are not extracting phase values from the speech envelope before downsampling.  
% Instead, we are downsampling the raw speech envelope and then analyzing the CSD between the EEG and the resampled speech envelope directly.
%   This should be fine, since CSD captures both shared power and phase relationships between the two signals at each frequency.  
% We don't need to extract phase values from the speech envelope in order to analyze its relationship with the EEG in the frequency domain.

% unlike our phase coherence analyis, we are not band-pass filtering the speech envelope before downsampling.
% because CSD captures both shared power and phase relationships between the two signals at each frequency, 
% we can analyze the relationship between the EEG and the raw speech envelope without band-pass filtering.

speech_envelope_resampled = resample(speech_envelope, num_samples_eeg, num_samples_speech);     % Downsample to match EEG sampling frequency
% The resampling function should produce a speech envelope that's the same length as the EEG.
% But in case it's off by a little bit, due to edge effects or the way the resampling function handles the input signal, 
% we should check the length of the resampled speech envelope and make sure it matches the number of samples in the EEG signal.
% we'll truncate or pad the resampled speech envelope as needed to ensure it has the same number of samples as the EEG signal, so that we can compute CSD without any issues.
if length(speech_envelope_resampled) > num_samples_eeg
    speech_envelope_resampled = speech_envelope_resampled(1:num_samples_eeg);  % truncate the speech if it's too long
elseif length(speech_envelope_resampled) < num_samples_eeg
    % If the resampled speech envelope is shorter than the EEG signal, we can pad it with zeros at the end to match the length.
    speech_envelope_resampled = [speech_envelope_resampled; zeros(num_samples_eeg - length(speech_envelope_resampled), 1)];
end

% let's plot the two signals to verify that they look like we expect.  
% The speech envelope should be a smoothed version of the original speech signal, and the EEG should contain both frequency components with some noise.
figure; 
subplot(3,1,1);
plot(t_speech, speech, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Speech Signal');
subplot(3,1,2);
plot(t_eeg, eeg, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Simulated EEG Signal');

% let's check the resampled speech envelope as well, to make sure it looks reasonable after downsampling.
subplot(3,1,3);
plot(t_eeg, speech_envelope_resampled, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Resampled Speech Envelope');

% let's check the length of the resampled speech envelope to make sure it matches the EEG length.
length(speech_envelope_resampled)

[CSD_one_sided, nfft] = CrossSpectralDensity(eeg, speech_envelope_resampled, num_samples_eeg, fs_eeg);
%%[CSD_one_sided, nfft] = CrossSpectralDensity(eeg, speech2_envelope, num_samples_eeg, fs_eeg);

% remember that the Fourier transform (and therefore CSD) of a real-valued signal is symmetric around the Nyquist frequency, 
% so we only need to analyze the first nfft/2 + 1 bins to get the full spectrum up to fs/2.
n_one_sided = nfft/2 + 1;  % number of frequency bins in the one-sided spectrum (including DC and Nyquist)
% Frequency axis for the one-sided spectrum
freqs = (0:n_one_sided-1) * (fs_eeg / nfft);
%% --- Extract magnitude and phase ---
CSD_mag = abs(CSD_one_sided);         % shared power at each frequency
CSD_phase = angle(CSD_one_sided);     % phase offset between x and y (radians)

%% --- Plot ---
figure;
subplot(2,1,1);
stem(freqs(1:100), CSD_mag(1:100), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|S_{xy}(f)|');
title('Cross-Spectral Density — Magnitude');
xlim([0 40]);
grid on;

% subplot(2,1,2);
% plot(freqs, CSD_phase * (180/pi), 'r', 'LineWidth', 1.5);
% xlabel('Frequency (Hz)');
% ylabel('Phase (degrees)');
% title('Cross-Spectral Density — Phase');
% xlim([0 fs/2]);
% ylim([-180 180]);
% grid on;

% let's plot phase in radians as well, since that's the native output of the CSD function and may be more intuitive for some analyses.
subplot(2,1,2);
plot(freqs, CSD_phase, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Cross-Spectral Density — Phase');
xlim([0 fs_eeg/2]);
ylim([-pi pi]);
grid on;

% given the signals I created above, why do I see a strong peak at 27 Hz?
% I think this is because the resampling process is introducing some artifacts into the speech envelope, 
% which are creating spurious frequency components that show up in the CSD.
% To test this, I can try analyzing the CSD between the EEG signal and the original speech signal (before resampling) to see if the 27 Hz peak is present there as well.
% If the 27 Hz peak is only present in the CSD with the resampled speech envelope and not in the CSD with the original speech signal, that would suggest that the resampling process is introducing artifacts that are creating this spurious peak.
% Alternatively, I can try analyzing the CSD between the EEG signal and a version of the speech envelope 
% that has been downsampled using a different method (e.g., simple decimation instead of resampling) to see if the 27 Hz peak is still present.


% the two signals I'm comparing above contain energy at 3 and 25 Hz, so I would expect to see peaks in the CSD magnitude at those frequencies.
% but i am also seeing a strong peak at 22 Hz, which is unexpected.  I need to investigate why that is happening.
% One possibility is that the resampling process is introducing some artifacts into the speech envelope, 
% which are creating spurious frequency components that show up in the CSD.
% how could that happen?
% The resampling process involves applying an anti-aliasing filter to the speech envelope before downsampling, 
% and if that filter is not designed properly or if there are edge effects, it could introduce artifacts into the resampled signal.

% can you suggest other possible reasons for the unexpected 22 Hz peak in the CSD between the EEG signal and the resampled speech envelope, besides artifacts from the resampling process?
% Other possible reasons for the unexpected 22 Hz peak in the CSD could include:
% 1) Spectral leakage: If the original speech signal contains strong frequency components that are not perfectly aligned with the FFT bins, 
% it can cause spectral leakage, which can create spurious peaks in the CSD. 
% This could be mitigated by using a longer window for the FFT or by applying a windowing function to the signals before computing the CSD.
% 2) Nonlinear interactions: If there are nonlinear interactions between the EEG signal and the speech envelope,
% it could create new frequency components that are not present in either signal alone, which could show up as unexpected peaks in the CSD.
% 3) Noise: If there is noise in either the EEG signal or the speech envelope, it could create spurious peaks in the CSD.
% To investigate the source of the unexpected 22 Hz peak, I could try analyzing the CSD between the EEG signal and 
% a version of the speech envelope that has been downsampled using a different method (e.g., simple decimation instead of resampling)
%  to see if the 22 Hz peak is still present.
%

