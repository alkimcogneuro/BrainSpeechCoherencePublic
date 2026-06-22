function [MSC_band, band_freqs, band_bin_idx] = BandAverageMSC(CoherenceResults, band_low, band_high)
    % =========================================================================================
    % Band-Averaged Magnitude Squared Coherence (MSC)
    % =========================================================================================
    % Collapse a CoherenceResults structure's per-frequency MSC_chanmeans down to a single MSC
    % value per channel, within a specified frequency band (e.g., 4-8 Hz).
    %
    % Arguments:
    %   CoherenceResults: a results structure produced by Apply2Dataset_CrossSpectralDensity,
    %                      containing (at minimum) the fields CSD_chanmeans, PSD_speech_chanmeans,
    %                      PSD_eeg_chanmeans (each [Num_channels x num_freqs]), and Freqs
    %                      ([1 x num_freqs] or [num_freqs x 1]).
    %   band_low:  lower edge of the frequency band of interest, in Hz (inclusive).
    %   band_high: upper edge of the frequency band of interest, in Hz (inclusive).
    %
    % Returns:
    %   MSC_band:     [Num_channels x 1] vector of band-averaged MSC values, one per channel.
    %   band_freqs:   the frequency bin values (Hz) that fell within [band_low, band_high] and
    %                 were used in the average. Returned for reference/sanity-checking.
    %   band_bin_idx: the indices into CoherenceResults.Freqs that were used. Returned for
    %                 reference/sanity-checking.
    %
    % -----------------------------------------------------------------------------------------
    % Design rationale: average spectra first, then take one ratio -- not average of ratios.
    % -----------------------------------------------------------------------------------------
    % MSC is a bounded ratio (0 to 1), not a power-like quantity, so the "correct" way to collapse
    % it across frequency bins is NOT to average MSC values directly. Instead, consistent with how
    % multi-trial averaging already works in this pipeline (CSD and PSD are averaged across trials
    % BEFORE dividing to get MSC_chanmeans), we average the complex CSD and the real-valued PSDs
    % across the frequency bins within the band first, and only then compute a single MSC value
    % from those band-averaged quantities:
    %
    %   CSD_band(ch)        = mean over band bins of CSD_chanmeans(ch, :)
    %   PSD_speech_band(ch) = mean over band bins of PSD_speech_chanmeans(ch, :)
    %   PSD_eeg_band(ch)    = mean over band bins of PSD_eeg_chanmeans(ch, :)
    %   MSC_band(ch)        = |CSD_band(ch)|^2 / (PSD_speech_band(ch) * PSD_eeg_band(ch))
    %
    % Averaging MSC ratios directly instead would treat every frequency bin as equally
    % informative regardless of how much power was actually present there, and would not be
    % consistent with the trial-averaging logic already used to compute MSC_chanmeans itself.
    % =========================================================================================

    % ---- Validate inputs -----------------------------------------------------------------
    required_fields = {'CSD_chanmeans', 'PSD_speech_chanmeans', 'PSD_eeg_chanmeans', 'Freqs'};
    for i = 1:numel(required_fields)
        if ~isfield(CoherenceResults, required_fields{i})
            error('BandAverageMSC:MissingField', ...
                'CoherenceResults is missing required field "%s". Was this produced by Apply2Dataset_CrossSpectralDensity?', ...
                required_fields{i});
        end
    end
    if band_low >= band_high
        error('BandAverageMSC:InvalidBand', ...
            'band_low (%.2f Hz) must be less than band_high (%.2f Hz).', band_low, band_high);
    end

    Freqs = CoherenceResults.Freqs(:)';  % ensure row vector, regardless of how it was stored

    % ---- Identify frequency bins within the requested band -------------------------------
    band_bin_idx = find(Freqs >= band_low & Freqs <= band_high);
    if isempty(band_bin_idx)
        error('BandAverageMSC:EmptyBand', ...
            ['No frequency bins found within [%.2f, %.2f] Hz. Available frequency range is ' ...
             '[%.2f, %.2f] Hz with %d bins. Check that the band falls within the analyzed range, ' ...
             'and that nfft/fs were set as expected when CoherenceResults was generated.'], ...
            band_low, band_high, Freqs(1), Freqs(end), numel(Freqs));
    end
    band_freqs = Freqs(band_bin_idx);

    % Warn (rather than silently proceeding) if the band is very coarsely sampled -- e.g., if
    % only one bin falls in the requested range, "band-averaging" reduces to that single bin,
    % which may not be what the user expects, especially for narrow bands or low-resolution nfft.
    if numel(band_bin_idx) == 1
        warning('BandAverageMSC:SingleBinBand', ...
            ['Only one frequency bin (%.3f Hz) falls within [%.2f, %.2f] Hz. The "band average" is ' ...
             'therefore just that single bin''s value. Consider whether nfft/frequency resolution ' ...
             'is fine enough for this band, or whether the band edges should be widened.'], ...
            band_freqs(1), band_low, band_high);
    end

    % ---- Average CSD and PSD across the band's bins, per channel -------------------------
    % CSD_chanmeans, PSD_speech_chanmeans, PSD_eeg_chanmeans are all [Num_channels x num_freqs].
    CSD_band        = mean(CoherenceResults.CSD_chanmeans(:, band_bin_idx), 2);         % complex, [Num_channels x 1]
    PSD_speech_band = mean(CoherenceResults.PSD_speech_chanmeans(:, band_bin_idx), 2);  % real,    [Num_channels x 1]
    PSD_eeg_band    = mean(CoherenceResults.PSD_eeg_chanmeans(:, band_bin_idx), 2);      % real,    [Num_channels x 1]

    % ---- Compute one MSC value per channel from the band-averaged spectra ----------------
    % MSC_band is a real-valued vector of length Num_channels, with one MSC value per channel, computed from the band-averaged CSD and PSD values for that channel.
    MSC_band = abs(CSD_band).^2 ./ (PSD_speech_band .* PSD_eeg_band);

    % ---- Flag (but do not silently hide) any NaN/Inf results -----------------------------
    % A NaN here means PSD_speech_band or PSD_eeg_band was zero (or both) for that channel,
    % which most likely indicates a degenerate trial somewhere upstream contaminated the
    % trial-averaged PSD/CSD for that channel. Surfacing this clearly here is much easier to
    % debug than discovering it later, downstream, in a group-level analysis.
    bad_channels = find(~isfinite(MSC_band));
    if ~isempty(bad_channels)
        warning('BandAverageMSC:NonFiniteResult', ...
            ['MSC_band is non-finite for %d channel(s): %s. This usually indicates zero or ' ...
             'near-zero band-averaged PSD for the speech or EEG signal in those channels, which ' ...
             'can happen if a degenerate trial (e.g., from a boundary or filtering issue) ' ...
             'contaminated the trial-averaged CSD/PSD for that channel.'], ...
            numel(bad_channels), mat2str(bad_channels(:)'));
    end
    %%call_topoplot_eeglab(MSC_band, CoherenceResults.Chanlocs, 'Mag Sq Coherence', 'MSC')
end
