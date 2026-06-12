function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, MSPC] = CrossSpectralDensity(sig_speech, sig_eeg, nfft)    
    % =========================================================
    % Cross-Spectral Density (CSD) Estimation
    % =========================================================
    % Calculate the cross spectral density between two signals (single trial). 
    % 
    % Function takes two arguments: 
    %   - sig_speech:  a vector containing the speech amplitude envelope samples
    %   - sig_eeg: a vector containing the EEG signal samples.
    %   - nfft:  the number of FFT points to use in the analysis.  This determines the frequency resolution of the CSD estimate.
    % 
    % Both signals should be sampled at the same rate and have the same length.
    % (any preprocessing steps to ensure this, such as resampling or segmenting, should be done before calling this function).
    % 
    % Returns:
    %  CSD_one_sided: the one-sided cross-spectral density estimate (complex-valued, length nfft/2 + 1)
    %  where nfft is the number of FFT points used in the analysis (equal to length fo the signal).
    % 
    % Mathematical Notes:
    % CSD of x and y is: S_xy(f) = X(f) * conj(Y(f)), 
    % where f indexes the frequency bins, X(f) and Y(f) are the Fourier Transforms of x and y, respectively at frequency f
    % 
    % At each frequency bin, CSD multiplies two complex values: 
    %       (A1 * e^(i theta1)) * (A2 * e^(i theta2))
    % which will multiply the magnitudes and add the phases.
    %       = A1*A2 * e^(i(theta1 + theta2))     
    % CSD is sensitive to both shared power and phase relationships between the two signals at each frequency.
    % CSD is a complex-valued function of frequency (one complex component per frequency bin):
    % - The magnitude |S_xy(f)| indicates how much power two signals x and y share at frequency f.
    % - The phase angle of S_xy(f) indicates the phase offset between the two signals at that frequency.
    % 
    % Why do we take the conjugate of Y(f) when computing the CSD?
    % conj(Y(f)) has a phase value that is the negation of the phase angle of Y(f). 
    % When we multiple X(f) by conj(Y(f)), we'll end up calculating the difference 
    % in phase between X and Y at each frequency, which is what we want to capture in the CSD. 
    % ----------------------------------------------------------------------------------------------
    % If signals are in row format, convert to column vector format.
    % Reason:  below, we apply the Hann window, and the window is a column vector.  
    % Element-wise multiplication of a row vector with a column vector would not work as intended.
    if isrow(sig_speech)
        sig_speech = sig_speech';
    end 
    if isrow(sig_eeg)
        sig_eeg = sig_eeg';
    end
    % Check that the two signals have the same length.
    signal_len = length(sig_speech);   % Calculate length of the signal (should be the same for both signals)
    if length(sig_eeg) ~= signal_len
        error('Input signals must have the same length. Length of sig_speech: %d, Length of sig_eeg: %d', signal_len, length(sig_eeg));
    end
    % ------------------------------------------------------------------------    
    % Create a Hann window of the same length as the signal, to taper the signal before computing the FFT.
    % Apply the Hann window to each signal by element-wise multiplication.
    % The windowing suppresses spectral leakage — the smearing of
    % energy from strong frequency components into neighboring bins.
    % ------------------------------------------------------------------------    
    win = hann(signal_len);  % Create a Hann window for our signal.
    sig_speech_win = win .* sig_speech;
    sig_eeg_win = win .* sig_eeg;
    % ------------------------------------------------------------------------    
    % Compute Discrete Fourier Transform (DFT) of each windowed segment using fft(). 
    % We'll set nfft (the number of FFT points) equal to signal length for simplicity.
    % ------------------------------------------------------------------------
    sig_speech_dft = fft(sig_speech_win, nfft);  
    sig_eeg_dft = fft(sig_eeg_win, nfft);
    % ------------------------------------------------------------------------
    % CSD of x and y is: S_xy(f) = X(f) * conj(Y(f)), 
    % The result of each multiplication is a complex value that captures both shared power and phase relationships 
    % between the two signals at each frequency.
    % ------------------------------------------------------------------------
    CSD = sig_speech_dft .* conj(sig_eeg_dft);
    % Normalize the CSD by the number of samples and the window power to get a proper estimate of the cross-spectral density.
    % The window power is the sum of the squared window values, 
    % which accounts for the energy reduction due to the tapering of the signal by the window.
    % 
    % Note that we normalize the values here, even though we will later normalize the CSD valuies by dividing by the geometric mean of the PSDs 
    % when we calculate the magnitude squared coherence (MSC).
    window_power = sum(win.^2); % Total power across the Hann window; use for normalization.
    CSD = CSD / (signal_len * window_power);    % Normalize

    % ------------------------------------------------------------------------
    % 
    % We'll calculate Power Spectral Density (PSD) as the CSD of each signal with itself.
    % ** Note ** when we multiply the DFT by its complex conjugate to calculate the PSD,
    % the result is a real-valued vector of power estimates at each frequency bin,
    % because the product of a complex number and its conjugate is always a real number 
    % with the imaginary phase values canceled out by the conjugate operation: 
    % i*theta - i*theta = 0.
    % ------------------------------------------------------------------------------
    PSD_speech = (sig_speech_dft .* conj(sig_speech_dft)) / (signal_len * window_power);  % PSD of the speech signal, normalized
    PSD_eeg = (sig_eeg_dft .* conj(sig_eeg_dft)) / (signal_len * window_power);  % PSD of the EEG signal, normalized.   
    % ---------------------------------------------------------------------------------
    % For CSD and PSD, we'll extract the one-sided spectra, up to the Nyquist.
    % For real-valued signals, we only need the first half of the spectrum 
    % We'll double all values, except for DC (bin 1) and Nyquist (bin end),
    % because the upper half of the spectrum is just the complex conjugate of the lower half, 
    % and we want to conserve total power in the one-sided spectrum.
    % ---------------------------------------------------------------------------------    
    n_one_sided = nfft/2 + 1;             % number of points in the one-sided spectrum (including DC and Nyquist)
    PSD_speech_one_sided = PSD_speech(1:n_one_sided);       % grab the one-sided PSD
    PSD_speech_one_sided(2:end-1) = 2 * PSD_speech_one_sided(2:end-1); % double values. 
    PSD_eeg_one_sided = PSD_eeg(1:n_one_sided);         % grab the one-side PSD
    PSD_eeg_one_sided(2:end-1) = 2 * PSD_eeg_one_sided(2:end-1); % double all values except for DC (bin 1) and Nyquist (bin end-1)
    CSD_one_sided = CSD(1:n_one_sided);
    CSD_one_sided(2:end-1) = 2 * CSD_one_sided(2:end-1); % double all values except for DC (bin 1) and Nyquist (bin end-1)    
    
    % Phase locking value (PLV).  
    % This is a measure of the consistency of the phase relationship between the two signals across trials.
    % turn each CSD vector into a unit vector by dividing by its magnitude.  
    % This removes information about magnitude, leaving only phase.
    % By averaging these unit vectors across trials, 
    % we can measure the consistency of the phase relationship between the two signals across trials, 
    % which is what PLV captures.
    %%% CSD_unit = CSD_one_sided ./ abs(CSD_one_sided); 
    % ---------------------------------------------------------------------------------
    % Calculate the SINGLE TRIAL magnitude squared PHASE coherence (MSPC) at each frequency. 
    % MSPC(f) = |S_xy(f)|^2 / (S_xx(f) * S_yy(f)), where S_xy is the CSD, and S_xx and S_yy are the PSDs of the two signals.
    % We're normalizing the magnitude squared CSD by the geometric mean of the two individual power spectra.
    % This is a value between 0 and 1, where 0 indicates no coherence (no consistent phase relationship) and 1 indicates perfect coherence (perfectly consistent phase relationship) at that frequency.
    % which is a measure of the consistency of the phase relationship between the two signals across trials.
    % I am going to call this "magnitude squared phase coherence", emphasizing the phase, 
    % because by normalizing the complex vector magnitude at the single trial level, 
    % we remove information about magnitude, leaving only phase.
    % So the coherence we measure by averaging these vectors is only about phase.  it's pure phase coherence.     %
    %
    % This measure is different from MSC, in which we will average the CSD values across trials and channels, 
    % and then calculate the MSC from those averages, which should give us an estimate of coherence
    % that reflects both consistency of phase and shared power across trials and channels.
    % ---------------------------------------------------------------------------------
    MSPC = abs(CSD_one_sided).^2 ./ (PSD_speech_one_sided .* PSD_eeg_one_sided);  % magnitude squared coherence at each frequency bin
end
