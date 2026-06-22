function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided] = CrossSpectralDensity(sig_speech, sig_eeg, nfft, fs)    
    % =========================================================
    % Cross-Spectral Density (CSD) Estimation
    % =========================================================
    % Calculate the cross spectral density between two signals (single trial). 
    % 
    % Function takes two arguments: 
    %   - sig_speech:  a vector containing the speech amplitude envelope samples
    %   - sig_eeg: a vector containing the EEG signal samples.
    %   - nfft: the number of FFT points to use in the analysis.  This determines the frequency resolution of the CSD estimate.
    %   - fs:  the sample rate of the signals.
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
    if nfft < signal_len
        warning('nfft (%d) is less than signal length (%d); signal will be truncated.', nfft, signal_len);
    elseif nfft > signal_len
        warning('nfft (%d) is greater than signal length (%d); signal will be zero-padded.', nfft, signal_len);
    end
    
    % mean-center the input signals
    % This will reduce a spike in energy for the DC offset component in the FFT. 
    sig_speech = sig_speech - mean(sig_speech);
    sig_eeg = sig_eeg - mean(sig_eeg);
    % ------------------------------------------------------------------------    
    % Create a Hann window of the same length as the signal, to taper the signal before computing the FFT.
    % Apply the Hann window to each signal by element-wise multiplication.
    % The windowing suppresses spectral leakage — the smearing of
    % energy from strong frequency components into neighboring bins.
    % ------------------------------------------------------------------------    
    win = hann(signal_len);  % Create a Hann window for our signal.
    window_power = sum(win.^2); % Total power across the Hann window; use for normalization.
                                % The window power is the sum of the squared window values, 
                                % which accounts for the energy reduction due to the tapering of the signal by the window.

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
    % "power per Hz"
    % 
    % Note that we normalize the values here, even though we will later normalize the CSD values by dividing by the geometric mean of the PSDs 
    % when we calculate the magnitude squared coherence (MSC). That's ok, since we'll normalize the numerator and denominator of the MSC by the same factor.      
    CSD = CSD / (fs * window_power);
    % ------------------------------------------------------------------------
    % We'll calculate Power Spectral Density (PSD) as the CSD of each signal with itself.
    % ** Note ** when we multiply the DFT by its complex conjugate to calculate the PSD,
    % the result is a real-valued vector of power estimates at each frequency bin,
    % because the product of a complex number and its conjugate is always a real number 
    % with the imaginary phase values canceled out by the conjugate operation: 
    % i*theta - i*theta = 0.
    % ------------------------------------------------------------------------------
    PSD_speech = (sig_speech_dft .* conj(sig_speech_dft)) / (fs * window_power);  % PSD of the speech signal, normalized
    PSD_eeg = (sig_eeg_dft .* conj(sig_eeg_dft)) / (fs * window_power);  % PSD of the EEG signal, normalized.   
    % ---------------------------------------------------------------------------------
    
    % ---------------------------------------------------------------------------------
    % For CSD and PSD, we'll extract the one-sided spectra, up to (and including) the
    % highest unique frequency bin.
    % For a real-valued signal of length nfft, the spectrum is conjugate-symmetric: bin k
    % and bin (nfft-k) are complex conjugates of one another. Bin 0 (DC) is always its own
    % mirror (real-valued) and is never doubled. Whether there is a second "self-mirrored"
    % bin depends on the parity of nfft:
    %   - nfft EVEN: bin nfft/2 (the Nyquist bin) maps to itself and is also never doubled.
    %     The one-sided spectrum has nfft/2 + 1 bins (indices 0 .. nfft/2); we double every
    %     bin strictly between DC and Nyquist.
    %   - nfft ODD: there is no exact Nyquist bin -- every bin from 1 up to (nfft-1)/2 has a
    %     distinct conjugate partner elsewhere in the spectrum. The one-sided spectrum has
    %     (nfft+1)/2 bins (indices 0 .. (nfft-1)/2); we double every bin EXCEPT DC, including
    %     the last one. (Using the even-case logic here would wrongly skip doubling the last
    %     bin, undercounting power at the highest analyzed frequency.)
    % In both cases we double all bins except DC and, only when nfft is even, the final
    % (Nyquist) bin, so that summing power across the one-sided spectrum still equals the
    % total power in the full two-sided spectrum.
    % These return objects are all [1 x n_one_sided]
    % ---------------------------------------------------------------------------------    
    if mod(nfft, 2) == 0                    % if nfft is even
        n_one_sided = nfft/2 + 1;          % bins 0 .. nfft/2 (DC through Nyquist), inclusive
        doubled_bins = 2:(n_one_sided - 1); % every bin except DC (1st) and Nyquist (last)
    else
        n_one_sided = (nfft + 1)/2;        % bins 0 .. (nfft-1)/2; no exact Nyquist bin exists
        doubled_bins = 2:n_one_sided;       % every bin except DC (1st); last bin IS doubled
    end
    PSD_speech_one_sided = PSD_speech(1:n_one_sided);       % grab the one-sided PSD
    PSD_speech_one_sided(doubled_bins) = 2 * PSD_speech_one_sided(doubled_bins);
    PSD_eeg_one_sided = PSD_eeg(1:n_one_sided);         % grab the one-side PSD
    PSD_eeg_one_sided(doubled_bins) = 2 * PSD_eeg_one_sided(doubled_bins);
    CSD_one_sided = CSD(1:n_one_sided);
    CSD_one_sided(doubled_bins) = 2 * CSD_one_sided(doubled_bins);
    
end
