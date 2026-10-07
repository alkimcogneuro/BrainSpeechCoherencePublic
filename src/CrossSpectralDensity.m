function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided] = CrossSpectralDensity(sig_speech, sig_eeg, nfft, fs)    
    % =========================================================
    % Cross-Spectral Density (CSD) Estimation
    % =========================================================
    % Calculate the cross spectral density between two signals at the level of a single trial. 
    % 
    % Function takes two arguments: 
    %   - sig_speech:  a vector containing the speech amplitude envelope samples
    %   - sig_eeg: a vector containing EEG signal samples.
    %   - nfft: the number of FFT points to use in the analysis.  This determines the frequency resolution of the CSD estimate.
    %   - fs:  the sample rate of the signals.
    % Both signals should be sampled at the same rate and have the same length.
    % 
    % Any preprocessing steps to ensure this, such as filtering, resampling or segmenting, 
    % should be done before calling this function.
    % 
    % Returns:
    %  CSD_one_sided: the one-sided cross-spectral density estimate (complex-valued, length nfft/2 + 1)
    %  where nfft is the number of FFT points used in the analysis (equal to length of the signal).
    % 
    % MATHEMATICAL NOTES:
    % CSD is a complex-valued function of frequency (one complex component per frequency bin):
    % CSD of x and y is: S_xy(f) = X(f) * conj(Y(f)), 
    % where f indexes the frequency bins, X(f) and Y(f) are the Fourier Transforms of x and y, respectively at frequency f
    % and conj(Y(f)) is the complex conjugate of Y(f).
    %
    % At each frequency bin, CSD multiplies two complex values: 
    %       (A1 * e^(i theta1)) * (A2 * e^(i theta2))
    % which will multiply the magnitudes and add the phases,
    % yielding a complex value at that frequency bin:
    %       = A1*A2 * e^(i(theta1 + theta2))     
    % 
    % At each frequency bin, the CSD value is sensitive to both shared power and phase relationships between the two signals.
    % - The magnitude |S_xy(f)| indicates how much power two signals x and y share at frequency f.
    % - The phase angle of S_xy(f) indicates the phase offset between the two signals at that frequency.
    % 
    % Conjugate flips the signs on the imaginary part:  
    % Why do we take the conjugate of Y(f) when computing the CSD?
    % conj(Y(f)) negates the phase values of Y(f). 
    % Multiplying X(f) and conj(Y(f)) adds the phase of X(f) to the negated phase of Y(f).
    % which calculate the difference in phase between the two signals at each frequency bin
    % which is what we want to capture in the CSD. 
    % ----------------------------------------------------------------------------------------------
    % 
    % ----------------------------------------
    % PREPROCESSING STEPS.
    % ----------------------------------------
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
    % Originally, the signals have different sampling rates, 
    % so we should resample the speech signal to match the EEG sampling rate before calling this function.
    signal_len = length(sig_speech);   % Calculate length of the signal (should be the same for both signals)
    if length(sig_eeg) ~= signal_len
        error('Input signals must have the same length. Length of sig_speech: %d, Length of sig_eeg: %d', signal_len, length(sig_eeg));
    end
    if nfft < signal_len
        warning('nfft (%d) is less than signal length (%d); signal will be truncated.', nfft, signal_len);
    elseif nfft > signal_len
        warning('nfft (%d) is greater than signal length (%d); signal will be zero-padded.', nfft, signal_len);
    end
    
    % I used to mean-center the signals, but that has now been replaced by 
    % a prior pre-processing step that z-scores the signals before calling this function.
    % Mean-center the input signals
    % This will reduce a spike in energy for the DC offset component in the FFT. 
    % sig_speech = sig_speech - mean(sig_speech);
    % sig_eeg = sig_eeg - mean(sig_eeg);
    % ------------------------------------------------------------------------    
    % Create a Hann window of the same length as the signals, to taper signals before computing the FFT.
    % Apply the Hann window to each signal by element-wise multiplication.
    % The windowing suppresses spectral leakage — the smearing of
    % energy from strong frequency components into neighboring bins.
    % ------------------------------------------------------------------------    
    win = hann(signal_len);  % Create a Hann window for our signal.
    window_power = sum(win.^2); % Total power across the Hann window; use for normalization.
                                % The window power is the sum of the squared window values, 
                                % which accounts for the energy reduction due to the tapering of the signal by the window.
    sig_speech_win = win .* sig_speech; % Apply Hann window
    sig_eeg_win = win .* sig_eeg;       % Apply Hann window
    % ------------------------------------------------------------------------    
    % Compute Discrete Fourier Transform (DFT) of each windowed segment using fft(). 
    % In most cases, we'll set nfft (the number of FFT points) equal to signal length for simplicity.
    % That will determine the value of nfft in the call to this function.
    % ------------------------------------------------------------------------
    sig_speech_dft = fft(sig_speech_win, nfft);  
    sig_eeg_dft = fft(sig_eeg_win, nfft);
    % ------------------------------------------------------------------------
    % Calculate CSD, as described above.
    % CSD of x and y is: S_xy(f) = X(f) * conj(Y(f)) 
    % ------------------------------------------------------------------------
    CSD = sig_speech_dft .* conj(sig_eeg_dft);
    % We will normalize the CSD by the number of samples and the window power to get a proper estimate of the cross-spectral density.
    % "power per Hz"
    % 
    % Note that we normalize the values here, even though we will later normalize the CSD values by dividing by the geometric mean of the PSDs 
    % when we calculate the magnitude squared coherence (MSC). That's ok, since we'll normalize the numerator and denominator of the MSC by the same factor.      
    CSD = CSD / (fs * window_power);
    % ------------------------------------------------------------------------
    % We'll calculate Power Spectral Density (PSD) of each signal as the CSD of each signal with itself.
    % Note: The PSD is a real-valued function of frequency:   
    % When we multiply the DFT by its complex conjugate, we add each phase value to its negated value, 
    % which cancels out -- i*theta - i*theta = 0.
    % This leaves only the squared magnitudes of the DFTs, at each frequency bin, which is a real number.
    % ------------------------------------------------------------------------------
    PSD_speech = (sig_speech_dft .* conj(sig_speech_dft)) / (fs * window_power);  % PSD of the speech signal, normalized
    PSD_eeg = (sig_eeg_dft .* conj(sig_eeg_dft)) / (fs * window_power);  % PSD of the EEG signal, normalized.   
    % ---------------------------------------------------------------------------------
    % ---------------------------------------------------------------------------------
    % ONE-SIDED SPECTRA
    % For CSD and PSD, we'll extract the one-sided spectra, 
    % up to and including the highest unique frequency bin.
    %
    % Math note:  
    % For a real-valued signal of length nfft, the spectrum is conjugate-symmetric: 
    % bin k and bin (nfft-k) are complex conjugates of one another, 
    % The right half of the spectrum is the (redundant) mirror image of the left half.
    % 
    % Because of the conjugate mirroring, we will keep only the left half, which contains all the unique frequency information, 
    % and we'll double the values of the bins that have a distinct conjugate partner in the right half of the spectrum.

    % Bin 0 (DC) is always its own mirror (real-valued). 
    % Whether there is a second "self-mirrored" bin depends on the parity of nfft:
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
    % The return objects are all [1 x n_one_sided]
    % ---------------------------------------------------------------------------------    
    if mod(nfft, 2) == 0                    % if nfft is even
        n_one_sided = nfft/2 + 1;           % bins 0 .. nfft/2 (DC through Nyquist), inclusive
        doubled_bins = 2:(n_one_sided - 1); % every bin except DC (1st) and Nyquist (last)
    else                                    % if nfft is odd.
        n_one_sided = (nfft + 1)/2;         % bins 0 .. (nfft-1)/2; no exact Nyquist bin exists
        doubled_bins = 2:n_one_sided;       % every bin except DC (1st); last bin IS doubled
    end
    PSD_speech_one_sided = PSD_speech(1:n_one_sided);       % grab the one-sided PSD
    PSD_speech_one_sided(doubled_bins) = 2 * PSD_speech_one_sided(doubled_bins);    % double the selected bins
    PSD_eeg_one_sided = PSD_eeg(1:n_one_sided);             % grab the one-side PSD
    PSD_eeg_one_sided(doubled_bins) = 2 * PSD_eeg_one_sided(doubled_bins);      % double the selected bins.
    CSD_one_sided = CSD(1:n_one_sided);
    CSD_one_sided(doubled_bins) = 2 * CSD_one_sided(doubled_bins);          % double the selected bins.
end
