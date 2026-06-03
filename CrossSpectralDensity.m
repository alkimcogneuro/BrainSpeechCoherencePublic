function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, MSPC] = CrossSpectralDensity(sig_speech, sig_eeg, nfft)    
    % =========================================================
    % Cross-Spectral Density (CSD) Estimation
    % =========================================================
    % Calculate the cross spectral density between two signals. 
    % This CSD calculation is for just one trial.
    % If we want to do multi-trial CSD averages, we'll average outputs from multiple calls to this function
    % 
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
    
    % Check that the two signals have the same length.
    % One potential reason for different lengths might be the resampling we do on the speech signal to match the EEG sampling rate.  
    % If there was an error in that resampling step, it could lead to a mismatch in lengths.
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
    % DFT produces one coefficient per frequency bin, 
    % where the bins correspond to frequencies from 0 to Nyquist frequency (inclusive) for the one-sided spectrum.
    % Each DFT coefficient is complex: magnitude = amplitude; angle = phase.
    % That is, each of coeffient is a vector in the complex plane, 
    % with a magnitude (distance from origin) and an angle (phase).
    % 
    % We'll set nfft (the number of FFT points) equal to signal length for simplicity.
    % If nfft > signal_len, the FFT will be zero-padded by Matlab's fft() function, 
    % which can make the spectrum look smoother but doesn't increase the true frequency resolution.
    % If nfft < signal_len, fft() would shorten the input signal to nfft samples, 
    % effectively cutting off the end of the segment, which would be bad.    
    % ------------------------------------------------------------------------
    sig_speech_dft = fft(sig_speech_win, nfft);  
    sig_eeg_dft = fft(sig_eeg_win, nfft);
    % ------------------------------------------------------------------------
    % Compute the CSD for the two signals. 
    % CSD of x and y is computed as S_xy(f) = X(f) * conj(Y(f)), 
    % where f indexes the frequency bins, X(f) and Y(f) are the DFT of x and y, respectively at frequency f
    % At each frequency bin, CSD multiplies two complex values: 
    %       (A1 * e^(i theta1)) * (A2 * e^(i theta2))
    % which will multiply the magnitudes and add the phases.
    %       = A1*A2 * e^(i(theta1 + theta2))     
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
    % Extract one-sided spectrum ---
    % For real-valued signals, we only need the first half of the spectrum produced by the DFT.
    % we call that the one-sided spectrum.
    % The full The DFT output is "two-sided" and symmetric for real signals.
    % The upper half of the spectrum (frequencies above Nyquist frequency) is just the complex conjugate of the lower half.
    % These are "mirror" frequencies that don't contain new information.
    % We keep only the positive frequencies (bins 1 to nfft/2 + 1)
    % and we **double the magnitude** to conserve total power.
    % but do not double the DC (bin 1) and Nyquist (bin nfft/2 + 1) components, which are unique and not mirrored.
    % The vector CSD_one_sided will be a complex-valued vector of length nfft/2 + 1, 
    % where each element corresponds to a specific frequency bin in the one-sided spectrum.  

    % what do we do here if nfft is an even number?  In that case, the Nyquist frequency is included in the one-sided spectrum and should not be doubled.  
    % The code below accounts for this by only doubling the values from bin 2 to bin nfft/2 (inclusive), 
    % and leaving bin 1 (DC) and bin nfft/2 + 1 (Nyquist) unchanged.
    % buf if nfft is even, then the index into CSD, nff2/2 + 1 will be a non-inter value, which will cause an error. 
    %  So we need to make sure nfft is an odd number, so that nfft/2 + 1 is an integer index into CSD.
    % what is normally done here?
    % In practice, nfft is often chosen to be a power of 2 for computational efficiency, which means it is often an even number.
    % In Matlab, when nfft is even, the one-sided spectrum includes the Nyquist frequency at index nfft/2 + 1, and that component should not be doubled.

    % we will also store the power spectral density (PSD) for the two signals
    % so that later, we can calculate the magnitude squared coherence (MSC) as MSC(f) = |S_xy(f)|^2 / (S_xx(f) * S_yy(f)), 
    % where S_xx and S_yy are the PSDs of x and y, respectively.
    % But we would need to be careful about how we calculate the PSDs, 
    % to make sure they are properly normalized and comparable to the CSD values.
    % In practice, the MSC is often calculated using Welch's method or other spectral estimation techniques that involve averaging across multiple segments of the data,
    % rather than using a single FFT of the entire signal, to get a more stable estimate of coherence.
    % but we are going to average the CSD values across trials and channels after we calculate them, 
    % and then we'll also average the PSD's across trials and channels, 
    % and then calculate the MSC from those averages, 
    % which should give us a stable estimate of coherence without needing to segment the data for Welch's method.
        
    % this seems like the wrong way to calculate the frequency vector, but it is actually correct
    % shouldn't we multiply the whole thing by the sampling rate to get frequencies in Hz?  No, because the sampling rate is not an argument to this function,
    %  and we want the function to be flexible and work with any sampling rate.  So we return the frequencies in cycles per sample, and then the caller can convert to Hz by multiplying by the sampling rate if they want.
    %  This will give us frequencies from 0 to Nyquist frequency (inclusive).
    % ------------------------------------------------------------------------

    n_one_sided = nfft/2 + 1;
    % ------------------------------------------------------------------------------
    % We'll calculate Power Spectral Density (PSD) 
    % The conventional way to calculate PSD would be to take 
    % the magnitude-squared of the DFT coefficients (the power), normalized in some appropriate way.
    % We'll calculate PSD by calculating the CSD of each signal with itself.
    % The CSD method is probably a little faster.
    % ** Note ** when we multiply the DFT by its complex conjugate to calculate the PSD,
    % the result is a real-valued vector of power estimates at each frequency bin,
    % because the product of a complex number and its conjugate is always a real number 
    % with the imaginary phase values canceled out by the conjugate operation: 
    % i*theta - i*theta = 0.
    %
    % We'll normalize the PSD signals by the length of the signal and the power of the window, 
    % just as we did for CSD.
    % 
    % The PSD values are POWER estimates at each frequency bin, 
    % which is what we need for calculating the magnitude squared coherence (MSC) later on.
    % when we use the geometric mean of the two PSDs in the denominator of the MSC calculation,
    % and we'll square the magnitude of the CSD in the numerator.
    % the CSD is also a power estimate, but it captures shared power between the two signals,
    % 
    % 
    % let's make sure this works the way we think it does.  
    % ------------------------------------------------------------------------------
    PSD_speech = (sig_speech_dft .* conj(sig_speech_dft)) / (signal_len * window_power);  % PSD of the speech signal, normalized
    PSD_eeg = (sig_eeg_dft .* conj(sig_eeg_dft)) / (signal_len * window_power);  % PSD of the EEG signal, normalized.   
    % ---------------------------------------------------------------------------------
    % For CSD and PSD, we'll grab just the one-sided spectra, up to the Nyquist 
    % We'll double all values, except for DC (bin 1) and Nyquist (bin end),
    % because the upper half of the spectrum is just the complex conjugate of the lower half, 
    % and we want to conserve total power in the one-sided spectrum.
    
    % ---------------------------------------------------------------------------------    
    PSD_speech_one_sided = PSD_speech(1:n_one_sided);       % grab the one-sided PSD
    PSD_speech_one_sided(2:end-1) = 2 * PSD_speech_one_sided(2:end-1); % double values. 
    PSD_eeg_one_sided = PSD_eeg(1:n_one_sided);         % grab the one-side PSD
    PSD_eeg_one_sided(2:end-1) = 2 * PSD_eeg_one_sided(2:end-1); % double all values except for DC (bin 1) and Nyquist (bin end-1)
    CSD_one_sided = CSD(1:n_one_sided);
    CSD_one_sided(2:end-1) = 2 * CSD_one_sided(2:end-1); % double all values except for DC (bin 1) and Nyquist (bin end-1)    
    
    % if i wanted to calculate power spectral density, shouldn't I take the magnitudes of values in PSD_speech and PSD_eeg?
    % Yes, the power spectral density (PSD) is typically calculated as the magnitude squared of the Fourier coefficients, 
    % which is what we have done here by multiplying the DFT by its complex conjugate.
    % But my vector PSD_speech_one_sided has not taken magnitudes by applying abs().  
    % It's only the product of the DFT and its conjugate, 
    % which should give us a real-valued vector of power values at each frequency bin,
    % because the product of a complex number and its conjugate is always a real number 
    % equal to the magnitude squared of the original complex number.
    % ahh, that's the part I didn't understand.  I thought that the product of the two complex numbers would be complex,
    % but it's real, because the phase values have been canceled out by the conjugate operation.  
    % So we do not need to take the magnitude of the PSD values, 
    % because they are already real-valued power estimates at each frequency bin.








    % Phase locking value (PLV).  
    % This is a measure of the consistency of the phase relationship between the two signals across trials.
    % turn each CSD vector into a unit vector by dividing by its magnitude.  
    % This removes information about magnitude, leaving only phase.
    % By averaging these unit vectors across trials, 
    % we can measure the consistency of the phase relationship between the two signals across trials, 
    % which is what PLV captures.
    CSD_unit = CSD_one_sided ./ abs(CSD_one_sided); 

    % ---------------------------------------------------------------------------------
    % Calculate the magnitude squared PHASE coherence (MSPC) at each frequency. 
    % MSPC(f) = |S_xy(f)|^2 / (S_xx(f) * S_yy(f)), where S_xy is the CSD, and S_xx and S_yy are the PSDs of the two signals.
    % So we're normalizing the magnitude squared CSD by the
    % geometric mean of the two individual power spectra.
    % The MSC is a value between 0 and 1, where 0 indicates no coherence (no consistent phase relationship) and 1 indicates perfect coherence (perfectly consistent phase relationship) at that frequency.
    % which is a measure of the consistency of the phase relationship between the two signals across trials.

    % ****************
    % ** question **:
    % 
    % What happens to the phases in the normalizing denominator term?
    % in the numerator, we've removed phase by grabbing the magnitude with abs()
    % But the PSDs are complex values, aren't they?
    % Is MSC made up of real or complex values?
    % I think we should extracting magnitudes from the PSDs, right?




    % I am going to call this "magnitude squared phase coherence", emphasizing the phase, 
    % because by normalizing the complex vector magnitude at the single trial level, 
    % we remove information about magnitude, leaving only phase.
    % So the coherence we measure by averaging these vectors is only about phase.  it's pure phase coherence. 
    %
    % We could call this Phase locking value (PLV).  maybe that's better than MSPC. 
    %
    % In a later calculation, we will average the CSD values across trials and channels, 
    % and then calculate the MSC from those averages, which should give us an estimate of coherence
    % that reflects both consistency of phase and shared power across trials and channels.
    % ---------------------------------------------------------------------------------
    
    MSPC = abs(CSD_one_sided).^2 ./ (PSD_speech_one_sided .* PSD_eeg_one_sided);  % magnitude squared coherence at each frequency bin

    % phase locking value.
    % Turn each CSD vector (single trial) into a unit vector by dividing by its magnitude


    %% If we want to extract magnitude and phase, we do this:
    % CSD_mag   = abs(CSD_one_sided);       % shared power at each frequency
    % CSD_phase = angle(CSD_one_sided);     % phase offset between x and y (radians)
end
