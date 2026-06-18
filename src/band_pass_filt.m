function signal_filtered_bandpass = band_pass_filt(signal, fs, highpass_cutoff, lowpass_cutoff)
    signal = signal(:);
    sig_len = length(signal);

    % Design filters using zero-pole-gain form, convert to SOS
    % This avoids the [b,a] numerical instability at low cutoff frequencies
    [z1, p1, k1] = butter(4, lowpass_cutoff/(fs/2), 'low');
    sos1 = zp2sos(z1, p1, k1);

    [z2, p2, k2] = butter(4, highpass_cutoff/(fs/2), 'high');
    sos2 = zp2sos(z2, p2, k2);

    % Estimate pad length based on the highpass cutoff (the slower-settling filter)
    % filtfilt recommends pad > 3x filter order; we use a frequency-aware estimate
    recommended_pad = round(5 * fs / highpass_cutoff);  % e.g. 10000 for 0.5 Hz
    pad_length = min(sig_len - 1, recommended_pad);

    % Warn if signal is too short for reliable filtering at this cutoff
    if pad_length < recommended_pad
        warning('band_pass_filt: signal may be too short for reliable filtering at %.2f Hz highpass. Consider longer epochs.', highpass_cutoff);
    end

    % Mirror pad
    signal_padded = [flip(signal(1:pad_length)); signal; flip(signal(end-pad_length+1:end))];

    % Apply filters using SOS (numerically stable)
    signal_filtered_lowpass  = filtfilt(sos1, 1, signal_padded);
    signal_filtered_bandpass = filtfilt(sos2, 1, signal_filtered_lowpass);

    % Trim padding
    signal_filtered_bandpass = signal_filtered_bandpass(pad_length+1:end-pad_length);
end