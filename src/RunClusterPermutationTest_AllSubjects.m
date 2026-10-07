function [ClusterResultsArray, SubjectSummary] = RunClusterPermutationTest_AllSubjects(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Run the Single-Subject Cluster-Based Permutation Test Across All Subjects
    % =========================================================================================
    % Thin batch wrapper around ClusterPermutationTest_ShuffleTest.m: runs it once per subject
    % in ShuffleTestResultsArray, using ONE shared channel-adjacency structure (built once, from
    % the first subject's Chanlocs, so every subject is tested against an identical montage
    % definition rather than silently rebuilding -- and potentially slightly re-deriving -- it
    % per subject), and rolls the results into a plain summary table you can eyeball or filter.
    %
    % This does NOT do any cross-subject statistical combination (no group-level test here) --
    % see open-items.md item 8 for why that's a separate, still-open design question when
    % clusters can sit at different channels per subject. This function only answers, per
    % subject, independently: does this subject show at least one significant cluster, and
    % where.
    %
    % Arguments:
    %   ShuffleTestResultsArray: struct array, one element per subject, as returned by
    %                            ShuffleTest_SpeechEEGCoherence.m. Build this the same way you
    %                            already do for GroupLevel_ShuffleTest_Analysis.m:
    %                                ShuffleTestResultsArray = [All_ShuffleTestResults.Coherence_Analysis_ShuffleTestResults];
    %
    % Optional name-value arguments (all passed straight through to
    % ClusterPermutationTest_ShuffleTest.m for every subject -- see that file for full detail):
    %   options.Threshold    (default norminv(0.95)): per-channel rank-Z cutoff for candidate
    %                        cluster membership.
    %   options.ClusterAlpha (default 0.05): cluster-mass significance threshold.
    %   options.Adjacency    (default []): reuse a precomputed adjacency matrix (e.g. one you
    %                        already built and checked with PlotChannelAdjacency.m) instead of
    %                        having this function build one from the first subject's Chanlocs.
    %                        Recommended once you've confirmed the auto-built one looks right --
    %                        build it once, look at it, then pass it in here explicitly.
    %
    % Returns:
    %   ClusterResultsArray: 1 x Num_subjects struct array; ClusterResultsArray(s) is exactly
    %                        what ClusterPermutationTest_ShuffleTest.m returns for subject s
    %                        (full per-cluster detail: channel membership, mass, p-value,
    %                        significance, the null distribution, etc.). Keep this around --
    %                        it's what a future topoplot-the-significant-clusters function will
    %                        need, not just the summary table below.
    %   SubjectSummary:      a table, one row per subject, with columns:
    %       Subj_id                 subject identifier (falls back to 'Subject_<n>' if a given
    %                                subject's ShuffleTestResults has no Subj_id field).
    %       NumClustersTotal        how many candidate clusters were found at all (significant
    %                                or not).
    %       NumSignificantClusters  how many of those cleared options.ClusterAlpha.
    %       HasSignificantCluster   logical -- true if NumSignificantClusters > 0. This is the
    %                                "does this subject show above-chance coherence somewhere"
    %                                flag.
    %       BestClusterMass, BestClusterPValue, BestClusterSize, BestClusterLabels:
    %                                detail for this subject's single largest-mass cluster
    %                                (whether or not it happens to be significant), so you can
    %                                see at a glance where/how strong the top candidate was even
    %                                for a subject with no significant cluster. BestClusterLabels
    %                                is a comma-separated string of that cluster's channel labels.
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Chanlocs consistency: this function requires every subject to have the same NUMBER of
    %    channels (same check GroupLevel_ShuffleTest_Analysis.m makes) and WARNS (does not error)
    %    if channel labels differ across subjects at the same index, since that would silently
    %    misalign "channel 12" between subjects if it ever happened, but might also just be
    %    harmless metadata variation -- worth knowing about either way, not necessarily fatal.
    % 2) Not run in an actual MATLAB session (none was reachable while writing this) --
    %    syntax-checked only, same caveat as ClusterPermutationTest_ShuffleTest.m itself.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.Threshold (1,1) double {mustBeReal} = norminv(0.95)
        options.ClusterAlpha (1,1) double {mustBeReal, mustBePositive} = 0.05
        options.Adjacency (:,:) logical = logical([])
    end

    Num_subjects = numel(ShuffleTestResultsArray);
    if Num_subjects < 1
        error('RunClusterPermutationTest_AllSubjects:NoSubjects', ...
            'ShuffleTestResultsArray is empty -- nothing to run.');
    end

    % ---- Chanlocs consistency checks (channel count required to match; labels just warned) ----
    Num_channels = numel(ShuffleTestResultsArray(1).True_MSC);
    Reference_labels = {ShuffleTestResultsArray(1).Chanlocs.labels};
    label_mismatch_subjects = {};
    for s = 2:Num_subjects
        if numel(ShuffleTestResultsArray(s).True_MSC) ~= Num_channels
            error('RunClusterPermutationTest_AllSubjects:ChannelCountMismatch', ...
                ['Subject %d has %d channels but subject 1 has %d. All subjects must share the ' ...
                 'same channel count/montage for a shared adjacency structure to make sense.'], ...
                s, numel(ShuffleTestResultsArray(s).True_MSC), Num_channels);
        end
        this_labels = {ShuffleTestResultsArray(s).Chanlocs.labels};
        if numel(this_labels) == numel(Reference_labels) && ~isequal(this_labels, Reference_labels)
            label_mismatch_subjects{end+1} = num2str(s); %#ok<AGROW>
        end
    end
    if ~isempty(label_mismatch_subjects)
        warning('RunClusterPermutationTest_AllSubjects:LabelMismatch', ...
            ['Channel labels at matching indices differ from subject 1 for subject(s): %s. ' ...
             'Channel counts match so this will still run, but double-check these subjects'' ' ...
             'channels are truly in the same order/montage as the others.'], ...
            strjoin(label_mismatch_subjects, ', '));
    end

    % ---- Shared adjacency structure (built once, reused for every subject) --------------------
    if isempty(options.Adjacency)
        fprintf('No adjacency supplied -- building one from subject 1''s Chanlocs.\n');
        Adjacency = BuildChannelAdjacency(ShuffleTestResultsArray(1).Chanlocs);
    else
        Adjacency = options.Adjacency;
        if ~isequal(size(Adjacency), [Num_channels, Num_channels])
            error('RunClusterPermutationTest_AllSubjects:AdjacencySizeMismatch', ...
                'options.Adjacency is %dx%d but there are %d channels.', ...
                size(Adjacency,1), size(Adjacency,2), Num_channels);
        end
    end

    % ---- Run the single-subject test for every subject -----------------------------------------
    Subj_id                = cell(Num_subjects, 1);
    NumClustersTotal        = zeros(Num_subjects, 1);
    NumSignificantClusters  = zeros(Num_subjects, 1);
    HasSignificantCluster   = false(Num_subjects, 1);
    BestClusterMass          = nan(Num_subjects, 1);
    BestClusterPValue        = nan(Num_subjects, 1);
    BestClusterSize          = nan(Num_subjects, 1);
    BestClusterLabels        = cell(Num_subjects, 1);

    for s = 1:Num_subjects
        if isfield(ShuffleTestResultsArray(s), 'Subj_id') && ~isempty(ShuffleTestResultsArray(s).Subj_id)
            this_id = ShuffleTestResultsArray(s).Subj_id;
        else
            this_id = sprintf('Subject_%d', s);
        end
        Subj_id{s} = this_id;

        fprintf('Running cluster permutation test for subject %d/%d (%s)...\n', s, Num_subjects, this_id);
        cr = ClusterPermutationTest_ShuffleTest(ShuffleTestResultsArray(s), ...
            'Threshold', options.Threshold, 'ClusterAlpha', options.ClusterAlpha, 'Adjacency', Adjacency);

        if s == 1
            ClusterResultsArray = cr;
        else
            ClusterResultsArray(s) = cr; %#ok<AGROW>
        end

        NumClustersTotal(s) = numel(cr.Clusters);
        NumSignificantClusters(s) = sum(cr.ClusterSignificant);
        HasSignificantCluster(s) = NumSignificantClusters(s) > 0;

        if NumClustersTotal(s) > 0
            [best_mass, best_idx] = max(cr.ClusterMass);
            BestClusterMass(s)   = best_mass;
            BestClusterPValue(s) = cr.ClusterPValue(best_idx);
            BestClusterSize(s)   = cr.ClusterSize(best_idx);
            BestClusterLabels{s} = strjoin(cr.ClusterLabels{best_idx}, ', ');
        else
            BestClusterLabels{s} = '';
        end

        if HasSignificantCluster(s)
            sig_idx = find(cr.ClusterSignificant);
            sig_descr = cell(1, numel(sig_idx));
            for k = 1:numel(sig_idx)
                sig_descr{k} = sprintf('[%s] mass=%.2f p=%.4f', ...
                    strjoin(cr.ClusterLabels{sig_idx(k)}, ','), cr.ClusterMass(sig_idx(k)), cr.ClusterPValue(sig_idx(k)));
            end
            fprintf('  -> SIGNIFICANT: %d cluster(s): %s\n', numel(sig_idx), strjoin(sig_descr, '; '));
        else
            fprintf('  -> no significant cluster (best candidate: mass=%.2f, p=%.4f, %d channel(s))\n', ...
                BestClusterMass(s), BestClusterPValue(s), BestClusterSize(s));
        end
    end

    SubjectSummary = table(Subj_id, NumClustersTotal, NumSignificantClusters, HasSignificantCluster, ...
        BestClusterMass, BestClusterPValue, BestClusterSize, BestClusterLabels);

    fprintf('\n%d/%d subjects have at least one significant cluster (alpha=%.3g).\n', ...
        sum(HasSignificantCluster), Num_subjects, options.ClusterAlpha);
end
