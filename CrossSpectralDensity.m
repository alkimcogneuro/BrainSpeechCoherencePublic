function [CSD_one_sided, CSD_freqs] = CrossSpectralDensity(sig_speech, sig_eeg, nfft)    
    % =========================================================
    % Cross-Spectral Density (CSD) Estimation
    % =========================================================
    % Calculate the cross spectral density between two signals. 
    % Function takes two arguments: 
    %   - sig_speech:  a vector containing the speech amplitude envelope samples
    %   - sig_eeg: a vector containing the EEG signal samples.
    % Both signals should be sampled at the same rate and have the same length.
    % (any preprocessing steps to ensure this, such as resampling or segmenting, should be done before calling this function).
    % 
    % Returns:
    %  CSD_one_sided: the one-sided cross-spectral density estimate (complex-valued, length nfft/2 + 1)
    %  where nfft is the number of FFT points used in the analysis (equal to length fo the signal).
    % 
    % Notes about cross spectral density (CSD):
    % CSD is sensitive to both shared power and phase relationships between the two signals at each frequency.
    % CSD is a complex-valued function of frequency (one complex component per frequency bin):
    % - The magnitude |S_xy(f)| indicates how much power two signals x and y share at frequency f.
    % - The phase angle of S_xy(f) indicates the phase offset between the two signals at that frequency.
    % CSD of x and y is computed as S_xy(f) = X(f) * conj(Y(f))
    % where X(f) and Y(f) are the Fourier transforms of x and y.
    % That is, we multiply X(f) by the complex conjugate of Y(f) at each frequency bin.
    % Put yet another way, the i-th coefficient of x's FFT is multiplied 
    % by the complex conjugate of the i-th coefficient of y's FFT.
    % 
    % NOTE ABOUT THE MATH:
    % Why do we take the conjugate of Y(f) when computing the CSD?
    % The conjugate on Y is what captures the phase difference between x and y at each frequency.
    % The conjugate operation negates the phase angle of Y(f), 
    % which allows the product X(f) * conj(Y(f)) to correctly represent 
    % the phase relationship between the two signals in the CSD.
    % Note:  the complex conjugate of a complex number z = a + bi is conj(z) = a - bi, 
    % which negates the imaginary part and thus negates the phase angle.
    % When we multiply X(f) by conj(Y(f)), the resulting complex number's magnitude 
    % reflects the shared power at that frequency, 
    % and its angle reflects the phase DIFFERENCE  between the two signals at that frequency.
    % (algebraically, this happens because when we multiply two complex exponentials,
    % their magnitudes multiply and their angles add.  
    % By taking the conjugate of Y(f), we effectively subtract Y's angle from X's angle,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Another way to understand CSD is that it measures the cross-correlation of x and y in the frequency domain, 
    % 
    % Check if sig_speech/sig_eeg is a row or column vector and convert to column if necessary.
    % Reason:  below, we apply the Hann window, and the window is a column vector.  Element-wise multiplication of a row vector with a column vector would not work as intended.
    
    if isrow(sig_speech)
        sig_speech = sig_speech';
    end 
    if isrow(sig_eeg)
        sig_eeg = sig_eeg';
    end

    signal_len = length(sig_speech);   % calculate length of the signal (should be the same for both signals)
    
    % We'll set nfft (the number of FFT points) equal to signal length for simplicity.
    % If nfft > signal_len, the FFT will be zero-padded, 
    % which can make the spectrum look smoother but doesn't actually increase the true frequency resolution.
    % If nfft < signal_len, Matlab's fft() function would shorten the input signal to nfft samples, 
    % effectively cutting off the end of the segment, which would be bad.    
    win = hann(signal_len);  % create a Hann window for our signal.
    
    % Apply the Hann window to each signal by element-wise multiplication.
    % This tapers the signal to zero at its edges.
    % This suppresses spectral leakage — the smearing of
    % energy from strong frequency components into
    % neighboring bins.
    sig_speech_win = win .* sig_speech;
    sig_eeg_win = win .* sig_eeg;
    
    % Compute Discrete Fourier Transform (DFT) of each windowed segment using fft().
    % DFT outputs are complex: magnitude = amplitude, angle = phase.
    % More specifically, fft outputs a vector of complex numbers, 
    % where each element corresponds to a specific frequency bin.
    % Each of these elements is a vector in the complex plane, 
    % with a magnitude (distance from origin) and an angle (phase).
    sig_speech_dft = fft(sig_speech_win, nfft);
    sig_eeg_dft = fft(sig_eeg_win, nfft);
    
    % Compute the CSD for this segment: 
    % CSD of x and y is computed as S_xy(f) = X(f) * conj(Y(f)), 
    % where f indexes the frequency bins, X(f) is the DFT of x at frequency f, and Y(f) is the DFT of y at frequency f.
    CSD = sig_speech_dft .* conj(sig_eeg_dft);    
    % Normalize the CSD by the number of samples and the window power to get a proper estimate of the cross-spectral density.
    % The window power is the sum of the squared window values, which accounts for the energy reduction due to windowing.
    window_power = sum(win.^2); % total power across the Hann window; use for normalization.
    CSD = CSD / (signal_len * window_power);

    %% --- Extract one-sided spectrum ---
    % For real-valued signals, we only need the first half of the spectrum produced by the DFT.
    % we call that the one-sided spectrum.
    % The full The DFT output is "two-sided" and symmetric for real signals.
    % The upper half of the spectrum (frequencies above Nyquist frequency) is just the complex conjugate of the lower half.
    % These are "mirror" frequencies that don't contain new information.
    % We keep only the positive frequencies (bins 1 to nfft/2 + 1)
    % and we **double the magnitude** to conserve total power.
    % but do not double the DC (bin 1) and Nyquist (bin nfft/2 + 1) components, which are unique and not mirrored.
    % The vector CSD_ons_sided will be a complex-valued vector of length nfft/2 + 1, 
    % where each element corresponds to a specific frequency bin in the one-sided spectrum.  
    n_one_sided = nfft/2 + 1;
    CSD_one_sided = CSD(1:n_one_sided);
    CSD_one_sided(2:end-1) = 2 * CSD_one_sided(2:end-1); % double all values except for DC (bin 1) and Nyquist (bin end-1)    
    %% If we want to extract magnitude and phase, we do this:
    % CSD_mag   = abs(CSD_one_sided);       % shared power at each frequency
    % CSD_phase = angle(CSD_one_sided);     % phase offset between x and y (radians)
end
