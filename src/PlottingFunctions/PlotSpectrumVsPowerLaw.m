function fig = PlotSpectrumVsPowerLaw(Result, options)
    % ===============================================================================================================
    % Plot the output of SpectrumVsPowerLaw.m (or SpeechEnvelopeSpectrumVsPowerLaw.m).
    % ===============================================================================================================
    % Top panel (log-log): the observed PSD, the power-law fit 1/f^beta (solid where it was fitted, dotted where it is
    % an interpolation across the excluded band), and -- if the Monte-Carlo null was run -- the range that the fitted
    % power law would produce by chance (pointwise 95% envelope of Gaussian power-law noise analysed with the same
    % Welch settings), shaded in a pale tint of the fit line's color because it belongs to that line. The tested band
    % (darker gray) and the guard margin that is left out of the fit (lighter gray) are shaded and named in the legend.
    %
    % Bottom panel: the same information as a residual, 10*log10(observed / fit), in dB. Dots mark frequency bins whose
    % residual exceeds the FAMILY-WISE (max-statistic) threshold, i.e. that are individually significant after
    % correcting for having looked at every frequency in the fit range. The residual is the plot to read when
    % deciding whether the excess is a localized PEAK (rises and falls inside the band) or a broad KNEE-shaped hump
    % (see the caveats in SpectrumVsPowerLaw.m).
    %
    % Options: Title (auto if empty), YLabel, FreqLimitsHz ([lo hi]; default = half the lower fit edge .. 1.5x the
    % upper fit edge, capped at Nyquist), Figure (existing figure handle to draw into).
    % Returns the figure handle.
    % ===============================================================================================================
    arguments
        Result (1,1) struct
        options.Title char = ''
        options.YLabel char = 'PSD (units^2/Hz)'
        options.FreqLimitsHz (1,2) double = [NaN NaN]
        options.Figure = []
    end

    % ---- palette: blue = observed, orange = model (the fit line AND its chance range, drawn as a translucent tint of the
    % same orange), grays = context only (tested band, guard margin), dark gray = thresholds/text ----
    c_obs  = [42 120 214] / 255;
    c_fit  = [235 104 52] / 255;
    c_ref  = [82 81 78] / 255;
    c_null = c_fit;                 % patch color; FaceAlpha below turns it into a pale orange
    null_alpha = 0.28;
    c_band = [0.80 0.80 0.80];
    c_guard = [0.92 0.92 0.92];
    c_grid = [0.90 0.90 0.90];
    c_txt  = [82 81 78] / 255;

    f = Result.Freqs;
    xl = options.FreqLimitsHz;
    if any(isnan(xl))
        xl = [Result.FitRangeHz(1)/2, min(Result.FitRangeHz(2)*1.5, Result.fs/2)];
    end
    vis = f > 0 & f >= xl(1) & f <= xl(2);
    has_null = ~isempty(Result.Null);

    if isempty(options.Figure)
        fig = figure('Color', 'w', 'Position', [100 100 780 740]);
    else
        fig = figure(options.Figure); clf(fig);
        set(fig, 'Color', 'w');
    end
    ax1 = axes('Parent', fig, 'Position', [0.10 0.39 0.86 0.52]);
    ax2 = axes('Parent', fig, 'Position', [0.10 0.09 0.86 0.21]);

    % ---- data needed for axis limits -------------------------------------------------------------------------------
    vals = [Result.PSD(vis); Result.PSD_fit(vis)];
    if has_null
        null_hi = Result.PSD_fit .* 10.^(Result.Null.Residual_dB_hi / 10);
        null_lo = Result.PSD_fit .* 10.^(Result.Null.Residual_dB_lo / 10);
        vals = [vals; null_hi(vis & Result.RangeMask); null_lo(vis & Result.RangeMask)];
    end
    vals = vals(isfinite(vals) & vals > 0);
    yl = [min(vals) / 1.6, max(vals) * 1.6];

    excl = Result.ExcludeHz;   band = Result.BandHz;

    % ============================== TOP: log-log spectrum ==============================================================
    hold(ax1, 'on');
    set(ax1, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', yl);
    h_guard = patch('Parent', ax1, 'XData', [excl(1) excl(2) excl(2) excl(1)], 'YData', [yl(1) yl(1) yl(2) yl(2)], ...
        'FaceColor', c_guard, 'EdgeColor', 'none');
    h_band = patch('Parent', ax1, 'XData', [band(1) band(2) band(2) band(1)], 'YData', [yl(1) yl(1) yl(2) yl(2)], ...
        'FaceColor', c_band, 'EdgeColor', 'none');
    h_null = [];
    if has_null
        r = find(Result.RangeMask);
        h_null = patch('Parent', ax1, 'XData', [f(r); flipud(f(r))], 'YData', [null_lo(r); flipud(null_hi(r))], ...
            'FaceColor', c_null, 'EdgeColor', 'none', 'FaceAlpha', null_alpha);
    end
    r_all = Result.RangeMask;
    % fitted power law: solid on the fitted flanks, dotted through the excluded interval
    y_solid = Result.PSD_fit;  y_solid(~Result.FitMask) = NaN;
    dotted = conv(double(Result.ExcludeMask), [1 1 1], 'same') > 0 & r_all;
    y_dot = Result.PSD_fit;  y_dot(~dotted) = NaN;
    plot(ax1, f, y_dot, ':', 'Color', c_fit, 'LineWidth', 2);
    h_fit = plot(ax1, f, y_solid, '-', 'Color', c_fit, 'LineWidth', 2);
    h_obs = plot(ax1, f(vis), Result.PSD(vis), '-', 'Color', c_obs, 'LineWidth', 2);

    text(sqrt(band(1) * band(2)), yl(2) / 1.25, sprintf('%g-%g Hz', band(1), band(2)), 'Parent', ax1, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', c_txt, 'FontSize', 9);

    excess_txt = sprintf('excess in band: %+.1f dB', Result.BandExcess_dB);
    if ~isempty(Result.Bootstrap)
        excess_txt = sprintf('%s [%+.1f, %+.1f]', excess_txt, Result.Bootstrap.CI_dB(1), Result.Bootstrap.CI_dB(2));
    end
    stat_lines = {excess_txt};                 % cell array => one text line per element (portable across MATLAB releases)
    if has_null
        stat_lines{end+1} = sprintf('p = %.3g (Monte-Carlo)', Result.Null.PValue);
    end
    text(0.98, 0.96, stat_lines, 'Parent', ax1, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', c_txt, 'FontSize', 9);

    ylabel(ax1, options.YLabel);
    if Result.ExponentFixed
        fit_label = sprintf('Power-law fit: 1/f^{\\beta}, \\beta = %.2f (fixed)', Result.Exponent);
    else
        fit_label = sprintf('Power-law fit: 1/f^{\\beta}, \\beta = %.2f (fitted to flanks)', Result.Exponent);
    end
    % Concatenate (rather than index-assign): the handles are of different graphics types (lines and patches).
    handles = [h_obs, h_fit, h_band, h_guard];
    labels = {'Observed spectrum', fit_label, sprintf('Tested band (%g-%g Hz)', band(1), band(2)), ...
              'Left out of fit (band +/- guard)'};
    if ~isempty(h_null)
        handles = [handles, h_null];
        labels{end+1} = sprintf('Range expected by chance from fit (%g%%)', 100*(1 - Result.Options.Alpha)); %#ok<AGROW>
    end
    legend(ax1, handles, labels, 'Location', 'southwest', 'Box', 'off', 'FontSize', 9);

    ttl = options.Title;
    if isempty(ttl)
        ttl = sprintf('Spectrum vs. power-law null: %g-%g Hz band', band(1), band(2));
    end
    title(ax1, ttl, 'FontWeight', 'normal', 'Color', c_txt);

    % ============================== BOTTOM: residual in dB ============================================================
    res = Result.Residual_dB;
    ext = res(vis & Result.RangeMask);
    ylo = min([ext; -1]);  yhi = max([ext; 1]);
    if has_null
        ylo = min([ylo; Result.Null.Residual_dB_lo(Result.RangeMask)]);
        yhi = max([yhi; Result.Null.Residual_dB_hi(Result.RangeMask); Result.Null.MaxThreshold_dB]);
    end
    pad = 0.08 * (yhi - ylo);
    yl2 = [ylo - pad, yhi + pad];
    hold(ax2, 'on');
    set(ax2, 'XScale', 'log', 'XLim', xl, 'YLim', yl2);
    patch('Parent', ax2, 'XData', [excl(1) excl(2) excl(2) excl(1)], 'YData', [yl2(1) yl2(1) yl2(2) yl2(2)], ...
        'FaceColor', c_guard, 'EdgeColor', 'none');
    patch('Parent', ax2, 'XData', [band(1) band(2) band(2) band(1)], 'YData', [yl2(1) yl2(1) yl2(2) yl2(2)], ...
        'FaceColor', c_band, 'EdgeColor', 'none');
    if has_null
        r = find(Result.RangeMask);
        patch('Parent', ax2, 'XData', [f(r); flipud(f(r))], ...
            'YData', [Result.Null.Residual_dB_lo(r); flipud(Result.Null.Residual_dB_hi(r))], ...
            'FaceColor', c_null, 'EdgeColor', 'none', 'FaceAlpha', null_alpha);
        plot(ax2, [f(r(1)) f(r(end))], Result.Null.MaxThreshold_dB * [1 1], '--', 'Color', c_ref, 'LineWidth', 1);
        text(f(r(end)), Result.Null.MaxThreshold_dB, ' family-wise threshold', 'Parent', ax2, 'Color', c_txt, ...
            'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
    plot(ax2, xl, [0 0], '-', 'Color', c_txt, 'LineWidth', 0.5);
    plot(ax2, f(vis), res(vis), '-', 'Color', c_obs, 'LineWidth', 2);
    if has_null && any(Result.Null.FreqSignificant)
        sg = Result.Null.FreqSignificant & vis;
        plot(ax2, f(sg), res(sg), 'o', 'MarkerSize', 5, 'MarkerFaceColor', c_obs, 'MarkerEdgeColor', 'w', 'LineWidth', 0.75);
    end
    ylabel(ax2, 'Excess over fit (dB)');
    xlabel(ax2, 'Frequency (Hz)');

    % ---- shared cosmetics ---------------------------------------------------------------------------------------
    xt = [0.25 0.5 1 2 3 4 6 8 12 16 24 32 48 64];
    xt = xt(xt >= xl(1) & xt <= xl(2));
    for ax = [ax1, ax2]
        set(ax, 'XTick', xt, 'Box', 'off', 'TickDir', 'out', 'XColor', c_txt, 'YColor', c_txt, 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
        local_try_set(ax, 'GridColor', c_grid);
        local_try_set(ax, 'GridAlpha', 1);
        local_try_set(ax, 'MinorGridLineStyle', 'none');
    end
    set(ax1, 'XTickLabel', []);
    set(ax2, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), xt, 'UniformOutput', false));   % plain numbers, not 10^x
    set(fig, 'CurrentAxes', ax1);
end

function local_try_set(h, prop, val)
    try
        set(h, prop, val);
    catch
        % property not supported by this graphics system; cosmetic only
    end
end
