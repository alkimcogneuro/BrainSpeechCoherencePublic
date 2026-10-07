function [MaxStatResultsArray, SubjectSummary] = RunMaxStatisticPermutationTest_AllSubjects(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Run the Single-Subject Max-Statistic Permutation Test (No Clustering) Across All Subjects
    % =========================================================================================
    % Batch wrapper around MaxStatisticPermutationTest_ShuffleTest.m: runs it once per subject in
    % ShuffleTestResultsArray and rolls the results into a plain summary table. Added as a
    % simpler, spatially-blind counterpart to RunClusterPermutationTest_AllSubjects.m -- this
    % touches no channel-adjacency/spatial code at all, so it's a useful independent cross-check
    % while ClusterPermutationTest_ShuffleTest.m's clustering logic is still being verified: if
    % this finds clearly significant channels for subjects where the cluster test finds nothing,
    % or the "best channel" here doesn't show up anywhere in the cluster test's candidate
    % clusters for the same subject, that points at the clustering/adjacency step specifically
    % rather than at the underlying shuffle-test statistics (which both tests start from
    % identically).
    %
    % Arguments:
    %   ShuffleTestResultsArray: struct array, one element per subject, as returned by
    %                            ShuffleTest_SpeechEEGCoherence.m. Build it the same way as for
    %                            the other group-level/batch functions in this repo:
    %                                ShuffleTestResultsArray = [All_ShuffleTestResults.Coherence_Analysis_ShuffleTestResults];
    %
    % Optional name-value arguments (passed straight through to
    % MaxStatisticPermutationTest_ShuffleTest.m for every subject):
    %   options.Metric (default 'MSC'): 'MSC' or 'RankZ' -- see that file for the tradeoff.
    %   options.Alpha  (default 0.05): significance threshold.
    %
    % Returns:
    %   MaxStatResultsArray: 1 x Num_subjects struct array; MaxStatResultsArray(s) is exactly
    %                        what MaxStatisticPermutationTest_ShuffleTest.m returns for subject s
    %                        (full per-channel detail, not just the summary below).
    %   SubjectSummary:      a table, one row per subject, with columns:
    %       Subj_id                  subject identifier (falls back to 'Subject_<n>' if missing).
    %       HasSignificantChannel    logical -- true if at least one channel beats the
    %                                max-corrected null. This is the direct answer to "is there
    %                                one channel for which true MSC is significantly greater than
    %                                the null distribution of maximum MSC across shuffles."
    %       NumSignificantChannels   how many channels individually clear the max-corrected
    %                                threshold (can be 0 even when BestChannelPValue is close to
    %                                significant; can be >1 if several channels are all strong).
    %       BestChannelLabel, BestChannelStat, BestChannelPValue:
    %                                detail on this subject's single best channel, whether or not
    %                                it happened to be significant.
    %
    % Not run in an actual MATLAB session (none was reachable while writing this) --
    % syntax-checked only, same caveat as the rest of this analysis family.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.Metric (1,:) char {mustBeMember(options.Metric, {'MSC', 'RankZ'})} = 'MSC'
        options.Alpha (1,1) double {mustBeReal, mustBePositive} = 0.05
    end

    Num_subjects = numel(ShuffleTestResultsArray);
    if Num_subjects < 1
        error('RunMaxStatisticPermutationTest_AllSubjects:NoSubjects', ...
            'ShuffleTestResultsArray is empty -- nothing to run.');
    end

    Subj_id                = cell(Num_subjects, 1);
    HasSignificantChannel  = false(Num_subjects, 1);
    NumSignificantChannels = zeros(Num_subjects, 1);
    BestChannelLabel       = cell(Num_subjects, 1);
    BestChannelStat         = nan(Num_subjects, 1);
    BestChannelPValue       = nan(Num_subjects, 1);

    for s = 1:Num_subjects
        if isfield(ShuffleTestResultsArray(s), 'Subj_id') && ~isempty(ShuffleTestResultsArray(s).Subj_id)
            this_id = ShuffleTestResultsArray(s).Subj_id;
        else
            this_id = sprintf('Subject_%d', s);
        end
        Subj_id{s} = this_id;

        fprintf('Running max-statistic permutation test for subject %d/%d (%s)...\n', s, Num_subjects, this_id);
        ms = MaxStatisticPermutationTest_ShuffleTest(ShuffleTestResultsArray(s), ...
            'Metric', options.Metric, 'Alpha', options.Alpha);

        if s == 1
            MaxStatResultsArray = ms;
        else
            MaxStatResultsArray(s) = ms; %#ok<AGROW>
        end

        HasSignificantChannel(s)  = ms.HasSignificantChannel;
        NumSignificantChannels(s) = sum(ms.ChannelSignificant);
        BestChannelLabel{s}       = ms.BestChannelLabel;
        BestChannelStat(s)        = ms.BestChannelStat;
        BestChannelPValue(s)      = ms.BestChannelPValue;

        if HasSignificantChannel(s)
            fprintf('  -> SIGNIFICANT: best channel %s, stat=%.4f, p=%.4f (%d/%d channels individually significant)\n', ...
                ms.BestChannelLabel, ms.BestChannelStat, ms.BestChannelPValue, ...
                NumSignificantChannels(s), numel(ms.ChannelStat));
        else
            fprintf('  -> not significant (best channel %s, stat=%.4f, p=%.4f)\n', ...
                ms.BestChannelLabel, ms.BestChannelStat, ms.BestChannelPValue);
        end
    end

    SubjectSummary = table(Subj_id, HasSignificantChannel, NumSignificantChannels, ...
        BestChannelLabel, BestChannelStat, BestChannelPValue);

    fprintf('\n%d/%d subjects have at least one channel significant by the max-statistic test (alpha=%.3g, metric=%s).\n', ...
        sum(HasSignificantChannel), Num_subjects, options.Alpha, options.Metric);
end
