function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, NumSegments] = WelchCrossSpectralDensity(sig_speech, sig_eeg, fs, window_len_samples, overlap_fraction)
    % =========================================================
    % Welch-style (segmented, overlapping-window) Cross-Spectral Density Estimation
    % for one speech envelope signal and one EEG channel signal.
    % =========================================================
    % Splits each input signal into overlapping windows of length=window_len_samples.
    % Computes CSD/PSD for each window using CrossSpectralDensity.m
    % (called once per segment, with nfft == window_len_samples, so each segment is windowed,
    % FFT'd, and normalized exactly the way CrossSpectralDensity.m already does for a
    % whole-trial estimate), and averages the resulting one-sided CSD/PSD across all segments.
    %
    % WHY use the Welch's windowed method?
    % A single Hann-tapered periodogram over an entire trial assumes one stable phase
    % relationship holds for the whole trial. Real speech-EEG coupling isn't guaranteed to hold
    % a single phase relationship for many seconds at a stretch, so instead of one long window
    % per trial, this computes several shorter-window estimates and averages them.
    % This is consistent with how ComputeCoherenceForPairing.m already averages CSD/PSD across TRIALS.
    % before taking the MSC ratio (never averages ratios). 
    % This function does the same thing one level down, across SEGMENTS within a trial: 
    % average the complex CSD and the real PSDs across segments
    % first, then (downstream, in ComputeCoherenceForPairing.m) average again across trials, and
    % only then take the MSC ratio once, at the very end. 
    % Note: Averaging per-segment MSC ratios would be an error 
    %
    % ARGUMENTS:
    %   sig_speech, sig_eeg   : equal-length vectors (one trial's speech envelope / EEG channel).
    %   fs                    : sampling rate (Hz), shared by both signals.
    %   window_len_samples    : length of each analysis window, in samples. Must be <= length of
    %                           the input signals. Sets the frequency resolution (fs/window_len_samples
    %                           Hz per bin) and the minimum reliably-resolvable frequency (roughly
    %                           3 cycles / (window_len_samples/fs) Hz) -- see RecommendWelchWindowLength.m,
    %                           which should be used to choose this value per dataset/trial-length
    %                           rather than hardcoding a single window length across datasets whose
    %                           trials may be very different durations.
    %   overlap_fraction      : fraction of window_len_samples by which consecutive windows overlap,
    %                           in [0, 0.5]. Capped at 0.5 (50%) deliberately: segments overlapping
    %                           more than 50% share the majority of their samples with their
    %                           neighbor (a Hann window concentrates most of its weight toward the
    %                           segment's center, so heavily-overlapping segments are looking at
    %                           nearly the same data, not independent-ish observations), which would
    %                           overstate how much real averaging is happening without actually
    %                           providing the bias/variance benefit that's the point of Welch
    %                           segmentation. 50% is the standard convention
    %                           (Welch, 1967; also the default in MATLAB's pwelch/mscohere).
    %
    % Returns:
    %   CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided : same shape/meaning as
    %       CrossSpectralDensity.m's outputs. 
    %       But note that 
    %       each spectrum is now the MEAN across all Welch-segments' own one-sided spectra (complex mean for
    %       CSD, real mean for the PSDs).
    %       and the one-sided spectra have lengths that match the Welch-window length (window_len_samples/2+1,
    %       not the original signal length), 
    % 
    %   NumSegments : the number of segments actually averaged together, returned so callers can
    %       log/sanity-check how much averaging happened.
    % =========================================================

    if nargin < 5
        overlap_fraction = 0.5;
    end
    if overlap_fraction < 0 || overlap_fraction > 0.5
        error('WelchCrossSpectralDensity:OverlapOutOfRange', ...
            ['overlap_fraction must be in [0, 0.5] (0%% to 50%% overlap). Got %.4f. Overlap above 50%% is ' ...
             'deliberately not supported here -- consecutive Hann-windowed segments would share more than ' ...
             'half their samples and stop behaving like separate observations in the average. See this ' ...
             'function''s header, or RecommendWelchWindowLength.m, for the reasoning.'], overlap_fraction);
    end

    % ----------------------------------------------------------------------------------------
    % Row/column handling, matching CrossSpectralDensity.m's convention (that function does its
    % own row->column conversion too, so this is redundant but harmless -- kept here so this
    % function's own length/indexing logic below doesn't have to special-case orientation).
    % ----------------------------------------------------------------------------------------
    if isrow(sig_speech)
        sig_speech = sig_speech';
    end
    if isrow(sig_eeg)
        sig_eeg = sig_eeg';
    end

    signal_len = length(sig_speech);
    if length(sig_eeg) ~= signal_len
        error('WelchCrossSpectralDensity:LengthMismatch', ...
            'Input signals must have the same length. Length of sig_speech: %d, Length of sig_eeg: %d', ...
            signal_len, length(sig_eeg));
    end
    if window_len_samples > signal_len
        error('WelchCrossSpectralDensity:WindowLongerThanSignal', ...
            ['window_len_samples (%d) is longer than the input signal (%d samples) -- cannot fit even one ' ...
             'analysis window. Use a shorter window, or (if the trial itself is too short for the frequency ' ...
             'range you need) see RecommendWelchWindowLength.m.'], window_len_samples, signal_len);
    end

    % ------------------------------------------------------------------------------------------
    % Determine segment start indices. Step size is derived from the (<=50%-capped) overlap
    % fraction; segments are dropped from the trial's tail if they don't fully fit -- standard
    % Welch behavior, we never pad or wrap to force in a partial extra segment.
    % ------------------------------------------------------------------------------------------
    step_samples = round(window_len_samples * (1 - overlap_fraction));
    step_samples = max(step_samples, 1);  % guard against a rounding-to-zero step at extreme inputs

    NumSegments = floor((signal_len - window_len_samples) / step_samples) + 1;
    if NumSegments < 1
        error('WelchCrossSpectralDensity:NoSegmentsFit', ...
            'No analysis windows fit: signal_len=%d, window_len_samples=%d, step_samples=%d.', ...
            signal_len, window_len_samples, step_samples);
    end
    if NumSegments < 3
        warning('WelchCrossSpectralDensity:FewSegments', ...
            ['Only %d segment(s) fit in this %d-sample signal with a %d-sample window (step %d samples). ' ...
             'With this few segments, this estimate gets little of the averaging benefit Welch segmentation ' ...
             'is meant to provide over a single whole-trial estimate -- consider a shorter window if the ' ...
             'frequency range you need allows it. See RecommendWelchWindowLength.m.'], ...
            NumSegments, signal_len, window_len_samples, step_samples);
    end

    segment_starts = 1 + (0:(NumSegments-1)) * step_samples;

    % ------------------------------------------------------------------------------------------
    % Accumulate CSD/PSD across segments. Each segment is windowed, FFT'd, and normalized by
    % CrossSpectralDensity.m exactly as it already does for a whole-trial estimate (nfft ==
    % window_len_samples == the segment's own length, so no zero-padding/truncation warnings
    % from that function fire here).
    % ------------------------------------------------------------------------------------------
    CSD_sum = [];
    PSD_speech_sum = [];
    PSD_eeg_sum = [];
    for seg_idx = 1:NumSegments
        idx_range = segment_starts(seg_idx) : (segment_starts(seg_idx) + window_len_samples - 1);
        [csd_seg, psd_speech_seg, psd_eeg_seg] = CrossSpectralDensity( ...
            sig_speech(idx_range), sig_eeg(idx_range), window_len_samples, fs);

        if isempty(CSD_sum)
            CSD_sum = csd_seg;
            PSD_speech_sum = psd_speech_seg;
            PSD_eeg_sum = psd_eeg_seg;
        else
            CSD_sum = CSD_sum + csd_seg;
            PSD_speech_sum = PSD_speech_sum + psd_speech_seg;
            PSD_eeg_sum = PSD_eeg_sum + psd_eeg_seg;
        end
    end

    % Average the complex CSD and the real PSDs across segments -- NOT a per-segment MSC ratio
    % averaged afterward. See header.
    CSD_one_sided = CSD_sum / NumSegments;
    PSD_speech_one_sided = PSD_speech_sum / NumSegments;
    PSD_eeg_one_sided = PSD_eeg_sum / NumSegments;
end
