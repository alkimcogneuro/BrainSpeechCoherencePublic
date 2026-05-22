function signal_filtered_bandpass = band_pass_filt(signal, fs, highpass_cutoff, lowpass_cutoff)
    % Force to column vector
    signal = signal(:); 
    sig_len = length(signal);    
    % Use a pad length that isn't longer than the signal itself
    % 100 samples is usually enough for these frequencies, 
    % but we ensure it never exceeds sig_len.
    pad_length = min(sig_len - 1, 100); 
    
    % Mirror pad
    % We'll pad the signal by mirroring the start and end of the signal, 
    % to mitigate edge artifacts from filtering.
    signal_padded = [flip(signal(1:pad_length)); signal; flip(signal(end-pad_length+1:end))];     
    
    % Comment about the filtering we perform here:
    % butter returns the filter coefficients for a 4th order Butterworth filter.
    % this is a zero-phase filter, so it doesn't introduce a phase shift, which is important for our analysis.
    % It is applied using Matlab's filtfilt function, 
    % which applies the filter forward and backward to avoid phase distortion.
    % (the backward pass cancels out the phase shift introduced by the forward pass, 
    % resulting in zero phase distortion).
    % The lowpass filter is applied first to remove high-frequency noise, 
    % and then the highpass filter is applied to remove low-frequency drift, 
    % resulting in a band-pass filtered signal.
    % The cutoff frequencies are normalized by the Nyquist frequency (fs/2) 
    % because the butter function expects normalized frequencies.
    % It is an IIR filter, which is more computationally efficient than FIR filters for the same order, 
    % but it can have a non-linear phase response. 
    % However, since we are using filtfilt, the overall phase response is linearized.

    [b1, a1] = butter(4, lowpass_cutoff/(fs/2), 'low');
    signal_filtered_lowpass = filtfilt(b1, a1, signal_padded);
    
    [b2, a2] = butter(4, highpass_cutoff/(fs/2), 'high');
    signal_filtered_bandpass = filtfilt(b2, a2, signal_filtered_lowpass);
    
    % Trim padding
    signal_filtered_bandpass = signal_filtered_bandpass(pad_length+1:end-pad_length);   
end