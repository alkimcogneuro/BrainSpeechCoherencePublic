function [GroupResults] = GroupLevel_ShuffleTest_Analysis(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Group-Level Analysis of ShuffleTest_SpeechEEGCoherence Results Across Subjects
    % =========================================================================================
    % Takes an array of ShuffleTestResults structs -- one per subject, as produced by
    % ShuffleTest_SpeechEEGCoherence.m -- and tests, PER CHANNEL, whether coherence is reliably
    % greater than the shuffle-derived null ACROSS THE GROUP of subjects.
    %
    % UPDATE: this now standardizes each subject's True_MSC using ComputeRankBasedZScore.m (a
    % rank/percentile-based transform of True_MSC's position within Null_MSC), NOT the raw
    % ShuffleTestResults.ZScores field ((True_MSC - Null_Mean)/Null_Std). See
    % ComputeRankBasedZScore.m's header for the full rationale; short version: MSC is non-negative
    % and its shuffle-derived null distribution is not guaranteed to be symmetric/Gaussian, so the
    % raw Z-score's negative tail is mechanically compressed relative to its positive tail (True_MSC
    % can't go below 0, no matter how large Null_Std is) -- meaning a raw Z of -3 and +3 are NOT
    % equally-extreme evidence, even though they look that way on a symmetric scale. The rank-based
    % transform sidesteps this by working entirely from True_MSC's percentile rank within Null_MSC,
    % which is valid regardless of the null's shape. No shuffles need to be re-run for this --
    % it's computed post-hoc from the True_MSC/Null_MSC that ShuffleTest_SpeechEEGCoherence.m
    % already returns.
    %
    % Arguments:
    %   ShuffleTestResultsArray: a struct array, one element per subject, where each element is a
    %                            ShuffleTestResults struct as returned by
    %                            ShuffleTest_SpeechEEGCoherence.m (must contain at minimum the
    %                            fields .True_MSC, .Null_MSC, .Chanlocs, .Subj_id). Build this from
    %                            script_AnalyzeCoherence_withPermTest.m's All_ShuffleTestResults via:
    %                                ShuffleTestResultsArray = [All_ShuffleTestResults.Coherence_Analysis_ShuffleTestResults];
    %
    % Optional name-value arguments:
    %   options.FDR_alpha (default 0.05): alpha level used for the Benjamini-Hochberg FDR
    %                                      correction across channels.
    %   options.Tail (default 'right'): 'right', 'left', or 'both' -- passed to both ttest and
    %                                    signrank. See Design rationale below for why 'right' is
    %                                    the default here (unlike GroupLevel_MSC_Analysis.m, which
    %                                    leaves this at MATLAB's two-sided default).
    %
    % Returns:
    %   GroupResults: a structure containing:
    %     .RankZ          [Num_subjects x Num_channels] each subject's rank-based z-score (from
    %                     ComputeRankBasedZScore.m applied to that subject's True_MSC/Null_MSC),
    %                     stacked. Row order matches ShuffleTestResultsArray.
    %     .SubjectIDs     {Num_subjects x 1} cell array of Subj_id, same order.
    %     .GroupMean_RankZ [1 x Num_channels] mean RankZ across subjects, per channel.
    %     .GroupSEM_RankZ  [1 x Num_channels] SEM of RankZ across subjects, per channel.
    %     .tstat          [1 x Num_channels] one-sample t-statistic, per channel, testing
    %                     GroupMean_RankZ against 0.
    %     .p_ttest        [1 x Num_channels] one-sample t-test p-value, per channel.
    %     .p_signrank     [1 x Num_channels] one-sample Wilcoxon signed-rank test p-value
    %                     (nonparametric companion to the t-test), per channel.
    %     .p_ttest_fdr    [1 x Num_channels] Benjamini-Hochberg FDR-corrected p_ttest.
    %     .p_signrank_fdr [1 x Num_channels] Benjamini-Hochberg FDR-corrected p_signrank.
    %     .sig_ttest_fdr, .sig_signrank_fdr   [1 x Num_channels] logical, true where the
    %                     corresponding FDR-corrected p-value < FDR_alpha.
    %     .Num_subjects, .Num_channels, .Chanlocs, .Tail, .FDR_alpha
    %     .band_low, .band_high, .NumShuffles_perSubject   passthrough from the subjects, IF present
    %                     and consistent (see consistency check below); otherwise omitted.
    %
    % -----------------------------------------------------------------------------------------
    % Design rationale
    % -----------------------------------------------------------------------------------------
    % 1) One-sample, not paired: GroupLevel_MSC_Analysis.m runs a PAIRED test, because each
    %    subject there contributes two raw values per channel (a real-condition MSC and a
    %    control-condition MSC) that need to be compared to each other. Here, each subject has
    %    already been reduced to ONE number per channel (RankZ, which already encodes "this
    %    subject's real coherence vs. this subject's own shuffle-derived null, as a percentile").
    %    The group-level question is therefore simpler: is that per-subject summary number
    %    reliably different from zero across the sample of subjects? That's a one-sample test
    %    (against 0), not a paired test.
    % 2) Why this is the "random-effects" approach, and why it matters: this treats each subject's
    %    RankZ as one independent observation of a subject-level effect, and asks whether that
    %    effect generalizes across the sampled subjects -- the standard approach in group
    %    neuroimaging analysis (summary-statistics / random-effects group analysis), because it
    %    accounts for BETWEEN-SUBJECT variability, not just within-subject (permutation) noise.
    %    Simply averaging each subject's summary value together (as an earlier version of this
    %    pipeline did, first with -log10(p), then with the raw Z-score) does not do this -- it has
    %    no notion of whether subjects agree with each other, so one subject with an extreme value
    %    can dominate a plain average even if the other subjects show nothing. A channel where
    %    every subject has a small, consistent positive RankZ will have a LARGE t-statistic
    %    (reliable across subjects), while a channel where one subject has a huge RankZ and the
    %    rest are near zero will have a much smaller t-statistic despite a similar or even larger
    %    mean -- which is exactly the distinction a plain average erases and a t-test is built to
    %    capture.
    % 3) Tail = 'right' by default: ShuffleTest_SpeechEEGCoherence.m's own test is one-sided (does
    %    coherence exceed the shuffle-derived null?), so RankZ is only meaningfully interpretable
    %    in the "more coherence than chance" direction. Testing the group-level mean RankZ with a
    %    one-sided ('right') test carries that same directional hypothesis through to the group
    %    level, rather than silently reverting to a two-sided test (as GroupLevel_MSC_Analysis.m's
    %    inline ttest() call currently does, despite that file's own header describing a one-sided
    %    hypothesis -- worth reconciling there too, but out of scope for this file). Pass
    %    options.Tail = 'both' if you want the two-sided version instead.
    % 4) Both a t-test and a Wilcoxon signed-rank test are reported, for the same reason
    %    GroupLevel_MSC_Analysis.m reports both: with modest numbers of subjects, a nonparametric
    %    companion to the parametric test is good practice, since normality of RankZ across
    %    subjects (as opposed to within one subject's null distribution, which RankZ's own
    %    construction already sidesteps) is still an assumption, not a guarantee.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.FDR_alpha (1,1) double {mustBeReal, mustBePositive} = 0.05
        options.Tail (1,:) char {mustBeMember(options.Tail, {'right', 'left', 'both'})} = 'right'
    end

    % ---- Validate inputs --------------------------------------------------------------------
    Num_subjects = numel(ShuffleTestResultsArray);
    if Num_subjects < 2
        error('GroupLevel_ShuffleTest_Analysis:TooFewSubjects', ...
            'Need at least 2 subjects for a group-level analysis; got %d.', Num_subjects);
    end

    required_fields = {'True_MSC', 'Null_MSC', 'Chanlocs', 'Subj_id'};
    for s = 1:Num_subjects
        for f = 1:numel(required_fields)
            if ~isfield(ShuffleTestResultsArray(s), required_fields{f})
                error('GroupLevel_ShuffleTest_Analysis:MissingField', ...
                    ['Element %d of ShuffleTestResultsArray is missing required field "%s". ' ...
                     'Was this produced by ShuffleTest_SpeechEEGCoherence.m?'], s, required_fields{f});
            end
        end
    end

    Num_channels = numel(ShuffleTestResultsArray(1).Chanlocs);
    if numel(ShuffleTestResultsArray(1).True_MSC) ~= Num_channels
        error('GroupLevel_ShuffleTest_Analysis:ChanlocsMismatch', ...
            'Subject 1 (%s) has %d True_MSC values but %d Chanlocs entries. These must match.', ...
            ShuffleTestResultsArray(1).Subj_id, numel(ShuffleTestResultsArray(1).True_MSC), Num_channels);
    end
    for s = 2:Num_subjects
        if numel(ShuffleTestResultsArray(s).True_MSC) ~= Num_channels
            error('GroupLevel_ShuffleTest_Analysis:InconsistentChannelCount', ...
                ['Subject %d (%s) has %d channels, but subject 1 (%s) has %d. All subjects must ' ...
                 'have the same number of channels, in the same order, to combine into a ' ...
                 'group-level analysis.'], s, ShuffleTestResultsArray(s).Subj_id, ...
                numel(ShuffleTestResultsArray(s).True_MSC), ShuffleTestResultsArray(1).Subj_id, Num_channels);
        end
        if size(ShuffleTestResultsArray(s).Null_MSC, 2) ~= Num_channels
            error('GroupLevel_ShuffleTest_Analysis:InconsistentChannelCount', ...
                ['Subject %d (%s) has Null_MSC with %d columns, but expected %d channels (to match ' ...
                 'subject 1''s Chanlocs). Was this subject''s ShuffleTestResults generated with a ' ...
                 'different channel montage?'], s, ShuffleTestResultsArray(s).Subj_id, ...
                size(ShuffleTestResultsArray(s).Null_MSC, 2), Num_channels);
        end
    end

    % ---- Compute each subject's rank-based z-score, and stack across subjects -----------------
    % See ComputeRankBasedZScore.m for the full rationale on why this (rather than the raw
    % ShuffleTestResults.ZScores field) is the statistically appropriate per-subject quantity to
    % feed into a group-level test.
    RankZ = nan(Num_subjects, Num_channels);
    SubjectIDs = cell(Num_subjects, 1);
    for s = 1:Num_subjects
        RankZ(s, :) = ComputeRankBasedZScore(ShuffleTestResultsArray(s).True_MSC, ShuffleTestResultsArray(s).Null_MSC);
        SubjectIDs{s} = ShuffleTestResultsArray(s).Subj_id;
    end

    % ---- Check for non-finite RankZ values ------------------------------------------------------
    % ComputeRankBasedZScore.m already warns per-subject if it produces a non-finite result (which
    % should be rare given its continuity correction); surfacing that again here, aggregated across
    % subjects, mirrors GroupLevel_MSC_Analysis.m's handling of non-finite inputs.
    nan_mask = ~isfinite(RankZ);
    if any(nan_mask(:))
        [bad_subj, bad_chan] = find(nan_mask);
        warning('GroupLevel_ShuffleTest_Analysis:NonFiniteInputs', ...
            ['%d (subject, channel) entries in RankZ are non-finite (e.g., subject %d, channel ' ...
             '%d). These will be excluded per-channel (via omitnan / a finite-value mask) from the ' ...
             'group means, SEMs, and significance tests for the affected channels.'], ...
            numel(bad_subj), bad_subj(1), bad_chan(1));
    end

    % ---- Group descriptive statistics (mean, SEM across subjects, per channel) ----------------
    GroupMean_RankZ = mean(RankZ, 1, 'omitnan');
    GroupSEM_RankZ  = std(RankZ, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(RankZ), 1));

    % ---- One-sample significance testing against 0, per channel -------------------------------
    tstat      = nan(1, Num_channels);
    p_ttest    = nan(1, Num_channels);
    p_signrank = nan(1, Num_channels);

    for ch_idx = 1:Num_channels
        z_col = RankZ(:, ch_idx);
        valid_mask = isfinite(z_col);

        if sum(valid_mask) < 2
            % Not enough valid (finite) observations to run a test for this channel; leave as
            % NaN rather than letting ttest/signrank error out or silently mislead.
            continue;
        end

        [~, p_ttest(ch_idx), ~, stats] = ttest(z_col(valid_mask), 0, 'Tail', options.Tail);
        tstat(ch_idx) = stats.tstat;

        p_signrank(ch_idx) = signrank(z_col(valid_mask), 0, 'tail', options.Tail);
    end

    % ---- Multiple comparisons correction (Benjamini-Hochberg FDR), across channels ------------
    p_ttest_fdr    = benjamini_hochberg_fdr(p_ttest);
    p_signrank_fdr = benjamini_hochberg_fdr(p_signrank);

    sig_ttest_fdr    = p_ttest_fdr < options.FDR_alpha;
    sig_signrank_fdr = p_signrank_fdr < options.FDR_alpha;

    % ---- Assemble results structure ------------------------------------------------------------
    GroupResults = struct();
    GroupResults.RankZ      = RankZ;
    GroupResults.SubjectIDs = SubjectIDs;

    GroupResults.GroupMean_RankZ = GroupMean_RankZ;
    GroupResults.GroupSEM_RankZ  = GroupSEM_RankZ;

    GroupResults.tstat      = tstat;
    GroupResults.p_ttest    = p_ttest;
    GroupResults.p_signrank = p_signrank;

    GroupResults.p_ttest_fdr    = p_ttest_fdr;
    GroupResults.p_signrank_fdr = p_signrank_fdr;
    GroupResults.sig_ttest_fdr    = sig_ttest_fdr;
    GroupResults.sig_signrank_fdr = sig_signrank_fdr;

    GroupResults.Num_subjects = Num_subjects;
    GroupResults.Num_channels = Num_channels;
    GroupResults.Chanlocs = ShuffleTestResultsArray(1).Chanlocs;
    GroupResults.Tail = options.Tail;
    GroupResults.FDR_alpha = options.FDR_alpha;

    % Passthrough band/shuffle metadata, only if present AND consistent across subjects -- these
    % are for labeling/bookkeeping, not used in any computation above, so we degrade gracefully
    % (omit rather than error) if they're missing or inconsistent.
    if all(arrayfun(@(s) isfield(s, 'band_low'), ShuffleTestResultsArray)) && ...
            all(arrayfun(@(s) isfield(s, 'band_high'), ShuffleTestResultsArray)) && ...
            all(arrayfun(@(s) isfield(s, 'NumShuffles'), ShuffleTestResultsArray))
        band_lows  = arrayfun(@(s) s.band_low, ShuffleTestResultsArray);
        band_highs = arrayfun(@(s) s.band_high, ShuffleTestResultsArray);
        if all(band_lows == band_lows(1)) && all(band_highs == band_highs(1))
            GroupResults.band_low  = band_lows(1);
            GroupResults.band_high = band_highs(1);
            GroupResults.NumShuffles_perSubject = arrayfun(@(s) s.NumShuffles, ShuffleTestResultsArray);
        else
            warning('GroupLevel_ShuffleTest_Analysis:InconsistentBand', ...
                ['Subjects do not all share the same band_low/band_high. Omitting .band_low/.band_high ' ...
                 'from GroupResults -- check that all subjects were run with the same frequency band ' ...
                 'before combining them into a group-level analysis.']);
        end
    end
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function p_fdr = benjamini_hochberg_fdr(p_values)
    % Benjamini-Hochberg FDR correction across a vector of p-values (one per channel).
    % Duplicated from GroupLevel_MSC_Analysis.m's identical local helper -- see that file's copy
    % for the note on this being a candidate for a shared utility function if this pattern gets
    % duplicated further. NaN p-values (e.g., from channels with insufficient valid data) are left
    % as NaN and excluded from the correction procedure, rather than being treated as p=0 or p=1.
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
