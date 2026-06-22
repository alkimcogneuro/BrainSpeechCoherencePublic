function [GroupResults] = GroupLevel_MSC_Analysis(CoherenceResults_real, CoherenceResults_control, band_low, band_high, options)
    % =========================================================================================
    % Group-Level Magnitude Squared Coherence (MSC) Analysis
    % =========================================================================================
    % Combine band-averaged MSC values across multiple subjects, separately for a real-data
    % condition and a (e.g., RandomizeOnsets) control condition, and test -- per channel --
    % whether real-condition MSC is reliably greater than control-condition MSC across subjects.
    %
    % Arguments:
    %   CoherenceResults_real:    cell array (or struct array) of CoherenceResults structures,
    %                              one per subject, from the real-data condition.
    %   CoherenceResults_control: cell array (or struct array) of CoherenceResults structures,
    %                              one per subject, from the control condition (e.g., produced
    %                              with RandomizeOnsets=true). Must be the SAME LENGTH as
    %                              CoherenceResults_real, and in the SAME SUBJECT ORDER -- this
    %                              function assumes element i of each array corresponds to the
    %                              same subject, since the significance test is a PAIRED test.
    %   band_low:  lower edge of the frequency band of interest, in Hz (inclusive).
    %   band_high: upper edge of the frequency band of interest, in Hz (inclusive).
    %
    % Optional name-value arguments:
    %   options.FDR_alpha (default 0.05): alpha level used for the Benjamini-Hochberg FDR
    %                                      correction across channels.
    %
    % Returns:
    %   GroupResults: a structure containing:
    %     .MSC_real          [Num_subjects x Num_channels] band-averaged MSC, real condition
    %     .MSC_control       [Num_subjects x Num_channels] band-averaged MSC, control condition
    %     .MSC_diff          [Num_subjects x Num_channels] per-subject paired difference (real - control)
    %     .GroupMean_real     [1 x Num_channels] mean across subjects, real condition
    %     .GroupMean_control  [1 x Num_channels] mean across subjects, control condition
    %     .GroupSEM_real      [1 x Num_channels] SEM across subjects, real condition
    %     .GroupSEM_control   [1 x Num_channels] SEM across subjects, control condition
    %     .GroupMean_diff     [1 x Num_channels] mean of the paired differences
    %     .GroupSEM_diff      [1 x Num_channels] SEM of the paired differences
    %     .p_ttest            [1 x Num_channels] paired t-test p-value, per channel
    %     .tstat              [1 x Num_channels] paired t-test t-statistic, per channel
    %     .p_wilcoxon         [1 x Num_channels] Wilcoxon signed-rank test p-value, per channel
    %     .p_ttest_fdr        [1 x Num_channels] FDR (Benjamini-Hochberg) corrected p-values, t-test
    %     .p_wilcoxon_fdr     [1 x Num_channels] FDR (Benjamini-Hochberg) corrected p-values, Wilcoxon
    %     .sig_ttest_fdr      [1 x Num_channels] logical, true where p_ttest_fdr < FDR_alpha
    %     .sig_wilcoxon_fdr   [1 x Num_channels] logical, true where p_wilcoxon_fdr < FDR_alpha
    %     .Num_subjects, .band_low, .band_high, .band_freqs, .Chanlocs (if available), .FDR_alpha
    %
    % -----------------------------------------------------------------------------------------
    % Design rationale: paired test, not pooled.
    % -----------------------------------------------------------------------------------------
    % Each subject contributes one real-condition value and one control-condition value per
    % channel, and these two values share all of that subject's idiosyncratic characteristics
    % (electrode impedance, head anatomy, alertness, noise level, etc.). A paired test (testing
    % whether the per-subject real-minus-control DIFFERENCE is reliably greater than zero across
    % subjects) cancels out these subject-level baseline differences and is more statistically
    % powerful than treating the two conditions as independent (unpaired/pooled) samples.
    % We report both a paired t-test (parametric) and a Wilcoxon signed-rank test (nonparametric),
    % since MSC is a bounded [0,1] quantity and may not satisfy normality assumptions well,
    % especially with modest numbers of subjects.
    % =========================================================================================
    arguments
        CoherenceResults_real
        CoherenceResults_control
        band_low (1,1) double {mustBeReal}
        band_high (1,1) double {mustBeReal}
        options.FDR_alpha (1,1) double {mustBeReal, mustBePositive} = 0.05
    end

    % ---- Validate subject-level pairing ----------------------------------------------------
    Num_subjects = numel(CoherenceResults_real);
    if numel(CoherenceResults_control) ~= Num_subjects
        error('GroupLevel_MSC_Analysis:MismatchedSubjectCount', ...
            ['CoherenceResults_real has %d subjects but CoherenceResults_control has %d subjects. ' ...
             'These must be the same length and in the same subject order, since the significance ' ...
             'test is PAIRED per subject.'], Num_subjects, numel(CoherenceResults_control));
    end
    if Num_subjects < 2
        error('GroupLevel_MSC_Analysis:TooFewSubjects', ...
            'Need at least 2 subjects for a group-level analysis; got %d.', Num_subjects);
    end

    % ---- Helper to access either a cell array or a struct array uniformly ------------------
    get_subj = @(arr, idx) get_subject_struct(arr, idx);

    % ---- Band-average MSC for the first subject, to determine Num_channels and band_freqs -
    [MSC_band_first, band_freqs, ~] = BandAverageMSC(get_subj(CoherenceResults_real, 1), band_low, band_high);
    Num_channels = numel(MSC_band_first);

    MSC_real    = zeros(Num_subjects, Num_channels);
    MSC_control = zeros(Num_subjects, Num_channels);

    MSC_real(1, :) = MSC_band_first(:)';

    % ---- Band-average MSC for the control condition, subject 1 (for consistency check) -----
    [MSC_band_ctrl_first, band_freqs_ctrl_first, ~] = BandAverageMSC(get_subj(CoherenceResults_control, 1), band_low, band_high);
    if numel(MSC_band_ctrl_first) ~= Num_channels
        error('GroupLevel_MSC_Analysis:ChannelCountMismatch', ...
            ['Subject 1: real condition has %d channels but control condition has %d channels. ' ...
             'All CoherenceResults structures must have a consistent number of channels.'], ...
            Num_channels, numel(MSC_band_ctrl_first));
    end
    % Confirm the real and control conditions actually used the same frequency bins within the
    % band. If nfft differed between conditions (e.g., real and control epochs of different
    % lengths), the two band_freqs vectors could differ, which would make the paired comparison
    % below compare MSC values computed from different underlying frequency resolutions without
    % any warning. We check this explicitly rather than silently keeping only the real
    % condition's band_freqs for the output.
    if ~isequal(size(band_freqs), size(band_freqs_ctrl_first)) || ...
            max(abs(band_freqs - band_freqs_ctrl_first)) > 1e-9
        error('GroupLevel_MSC_Analysis:FrequencyAxisMismatch', ...
            ['Subject 1: the frequency bins within [%.2f, %.2f] Hz differ between the real and ' ...
             'control conditions. This likely means nfft (or fs) differed between how the real ' ...
             'and control CoherenceResults were generated, which would make a paired comparison ' ...
             'invalid. Check that both conditions used the same epoch length and sampling rate.'], ...
            band_low, band_high);
    end
    MSC_control(1, :) = MSC_band_ctrl_first(:)';

    % ---- Loop over remaining subjects, band-averaging real and control conditions ----------
    for subj_idx = 2:Num_subjects
        [MSC_band_real_i, ~, ~] = BandAverageMSC(get_subj(CoherenceResults_real, subj_idx), band_low, band_high);
        [MSC_band_ctrl_i, ~, ~] = BandAverageMSC(get_subj(CoherenceResults_control, subj_idx), band_low, band_high);

        if numel(MSC_band_real_i) ~= Num_channels
            error('GroupLevel_MSC_Analysis:ChannelCountMismatch', ...
                ['Subject %d (real condition) has %d channels, but subject 1 had %d channels. ' ...
                 'All subjects'' CoherenceResults must have the same number of channels to combine ' ...
                 'into a group-level analysis.'], subj_idx, numel(MSC_band_real_i), Num_channels);
        end
        if numel(MSC_band_ctrl_i) ~= Num_channels
            error('GroupLevel_MSC_Analysis:ChannelCountMismatch', ...
                ['Subject %d (control condition) has %d channels, but subject 1 had %d channels. ' ...
                 'All subjects'' CoherenceResults must have the same number of channels to combine ' ...
                 'into a group-level analysis.'], subj_idx, numel(MSC_band_ctrl_i), Num_channels);
        end

        MSC_real(subj_idx, :)    = MSC_band_real_i(:)';
        MSC_control(subj_idx, :) = MSC_band_ctrl_i(:)';
    end

    % ---- Check for subjects/channels with non-finite band-averaged MSC ---------------------
    % BandAverageMSC already warns per-subject if it produces non-finite values; here we check
    % whether that has left any NaNs that would silently propagate into (and invalidate) the
    % group-level statistics below, since mean()/std()/ttest() on data containing NaN will, by
    % default, either error or quietly return NaN themselves depending on the function. We
    % surface this explicitly rather than letting it happen implicitly downstream.
    nan_mask = ~isfinite(MSC_real) | ~isfinite(MSC_control);
    if any(nan_mask(:))
        [bad_subj, bad_chan] = find(nan_mask);
        warning('GroupLevel_MSC_Analysis:NonFiniteInputs', ...
            ['%d (subject, channel) entries are non-finite across the real/control band-averaged ' ...
             'MSC values (e.g., subject %d, channel %d). These will propagate as NaN into the ' ...
             'group means, SEMs, and significance tests for the affected channels. Consider ' ...
             'excluding the affected subjects/channels, or investigate the upstream cause ' ...
             '(see BandAverageMSC warnings for the affected subject(s)).'], ...
            numel(bad_subj), bad_subj(1), bad_chan(1));
    end

    % ---- Per-subject paired difference -------------------------------------------------------
    MSC_diff = MSC_real - MSC_control;  % [Num_subjects x Num_channels]

    % ---- Group descriptive statistics (mean, SEM across subjects, per channel) -------------
    GroupMean_real    = mean(MSC_real, 1, 'omitnan');
    GroupMean_control = mean(MSC_control, 1, 'omitnan');
    GroupMean_diff     = mean(MSC_diff, 1, 'omitnan');

    GroupSEM_real    = std(MSC_real, 0, 1, 'omitnan')    ./ sqrt(sum(isfinite(MSC_real), 1));
    GroupSEM_control = std(MSC_control, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(MSC_control), 1));
    GroupSEM_diff     = std(MSC_diff, 0, 1, 'omitnan')    ./ sqrt(sum(isfinite(MSC_diff), 1));

    % ---- Paired significance testing, per channel -------------------------------------------
    p_ttest    = nan(1, Num_channels);
    tstat      = nan(1, Num_channels);
    p_wilcoxon = nan(1, Num_channels);

    for ch_idx = 1:Num_channels
        real_col    = MSC_real(:, ch_idx);
        control_col = MSC_control(:, ch_idx);
        valid_mask  = isfinite(real_col) & isfinite(control_col);

        if sum(valid_mask) < 2
            % Not enough valid (finite) paired observations to run a test for this channel;
            % leave as NaN rather than letting ttest/signrank error out or silently mislead.
            continue;
        end

        % Paired t-test: tests whether the mean of (real - control) differs from zero.
        [~, p_ttest(ch_idx), ~, stats] = ttest(real_col(valid_mask), control_col(valid_mask));
        tstat(ch_idx) = stats.tstat;

        % Wilcoxon signed-rank test: nonparametric paired test, does not assume the
        % differences are normally distributed -- more robust given MSC is bounded [0,1].
        p_wilcoxon(ch_idx) = signrank(real_col(valid_mask), control_col(valid_mask));
    end

    % ---- Multiple comparisons correction (Benjamini-Hochberg FDR), across channels ---------
    p_ttest_fdr    = benjamini_hochberg_fdr(p_ttest);
    p_wilcoxon_fdr = benjamini_hochberg_fdr(p_wilcoxon);

    sig_ttest_fdr    = p_ttest_fdr < options.FDR_alpha;
    sig_wilcoxon_fdr = p_wilcoxon_fdr < options.FDR_alpha;

    % ---- Assemble results structure ----------------------------------------------------------
    GroupResults = struct();
    GroupResults.MSC_real    = MSC_real;
    GroupResults.MSC_control = MSC_control;
    GroupResults.MSC_diff    = MSC_diff;

    GroupResults.GroupMean_real    = GroupMean_real;
    GroupResults.GroupMean_control = GroupMean_control;
    GroupResults.GroupSEM_real     = GroupSEM_real;
    GroupResults.GroupSEM_control  = GroupSEM_control;
    GroupResults.GroupMean_diff    = GroupMean_diff;
    GroupResults.GroupSEM_diff     = GroupSEM_diff;

    GroupResults.p_ttest    = p_ttest;
    GroupResults.tstat      = tstat;
    GroupResults.p_wilcoxon = p_wilcoxon;

    GroupResults.p_ttest_fdr    = p_ttest_fdr;
    GroupResults.p_wilcoxon_fdr = p_wilcoxon_fdr;
    GroupResults.sig_ttest_fdr    = sig_ttest_fdr;
    GroupResults.sig_wilcoxon_fdr = sig_wilcoxon_fdr;

    GroupResults.Num_subjects = Num_subjects;
    GroupResults.Num_channels = Num_channels;
    GroupResults.band_low  = band_low;
    GroupResults.band_high = band_high;
    GroupResults.band_freqs = band_freqs;
    GroupResults.FDR_alpha = options.FDR_alpha;

    % Carry channel locations through for topographic plotting, if available on subject 1's
    % real-condition CoherenceResults (this field is populated by Apply2Dataset_CrossSpectralDensity).
    subj1_real = get_subj(CoherenceResults_real, 1);
    if isfield(subj1_real, 'Chanlocs')
        GroupResults.Chanlocs = subj1_real.Chanlocs;
    end
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function s = get_subject_struct(arr, idx)
    % Allow CoherenceResults_real / CoherenceResults_control to be passed in as either a cell
    % array of structs or a struct array, since the user indicated their storage format isn't
    % yet finalized. This indexes into either consistently.
    if iscell(arr)
        s = arr{idx};
    else
        s = arr(idx);
    end
end

function p_fdr = benjamini_hochberg_fdr(p_values)
    % Benjamini-Hochberg FDR correction across a vector of p-values (one per channel).
    % NaN p-values (e.g., from channels with insufficient valid data) are left as NaN and
    % excluded from the correction procedure, rather than being treated as p=0 or p=1.
    p_fdr = nan(size(p_values));
    valid_mask = isfinite(p_values);
    p_valid = p_values(valid_mask);
    m = numel(p_valid);
    if m == 0
        return;
    end

    [p_sorted, sort_idx] = sort(p_valid);
    ranks = (1:m);
    p_adj_sorted = p_sorted .* m ./ ranks;

    % Enforce monotonicity (standard step-up procedure): each adjusted p-value cannot be
    % smaller than the one after it in sorted order.
    p_adj_sorted = fliplr(cummin(fliplr(p_adj_sorted)));
    p_adj_sorted = min(p_adj_sorted, 1);  % clip at 1

    p_adj = nan(1, m);
    p_adj(sort_idx) = p_adj_sorted;

    p_fdr(valid_mask) = p_adj;
end
