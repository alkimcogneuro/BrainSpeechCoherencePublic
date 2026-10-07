function [ClusterResults] = ClusterPermutationTest_ShuffleTest(ShuffleTestResults, options)
    % =========================================================================================
    % Single-Subject Cluster-Based Permutation Test for Shuffle-Test MSC Results
    % =========================================================================================
    % Addresses inter-subject topographic variability in speech-EEG coherence (different
    % subjects may show real above-chance coherence at different, spatially separate channels,
    % e.g. because of functional-anatomy or electrode-placement differences) WITHOUT the
    % double-dipping problem you'd get from just picking "whichever channels look significant"
    % and then testing that same selection for significance again (see Kriegeskorte et al.,
    % 2009, "Circular analysis in systems neuroscience", for why that inflates false positives).
    %
    % This implements the standard cluster-based permutation approach (Maris & Oostenveld, 2007)
    % at the level of a SINGLE subject's shuffle test, reusing the null distribution you already
    % have (ShuffleTestResults.Null_MSC, [NumShuffles x Num_channels]) instead of requiring any
    % new data or a data split:
    %   1) Convert True_MSC to a per-channel rank-based Z (ComputeRankBasedZScore.m -- see that
    %      file for why raw Z is NOT appropriate here, given MSC's floor at 0).
    %   2) Threshold that Z at options.Threshold; group the surviving channels into spatially
    %      contiguous clusters using a channel-adjacency structure (BuildChannelAdjacency.m).
    %   3) Sum the per-channel Z within each cluster -> that cluster's "cluster mass".
    %   4) Build a null distribution of the LARGEST cluster mass achievable by chance: for each
    %      of the NumShuffles null draws already stored in Null_MSC, treat that single draw as a
    %      surrogate "true" value, compute ITS rank-based Z against the other NumShuffles-1 draws
    %      (leave-one-out -- same ComputeRankBasedZScore.m function, just called on a shuffle
    %      draw instead of the real True_MSC), threshold and cluster it the same way, and record
    %      the largest cluster mass found (0 if no cluster survives threshold).
    %   5) Each REAL cluster's p-value is how often the null's max-cluster-mass matched or beat
    %      that real cluster's mass: p = (1 + count) / (NumShuffles + 1) (same Phipson & Smyth,
    %      2010 estimator used throughout this codebase, e.g. ShuffleTest_SpeechEEGCoherence.m).
    %
    % Comparing every real cluster against a null built from the MAXIMUM cluster mass per
    % iteration (rather than each cluster against its own local null) is what makes this
    % family-wise-error corrected across the whole scan for multiple comparisons -- exactly the
    % same logic as Maris & Oostenveld's cluster-mass permutation test. Because every observed
    % cluster is tested against that same null, this does NOT limit you to a single "winning"
    % cluster per subject: if a subject genuinely shows two separate above-chance regions (e.g.
    % left AND right auditory cortex), and both clusters' masses beat the max-null distribution,
    % BOTH come back significant. "Max" is only used to build the null -- it is not a restriction
    % on how many real clusters can end up significant.
    %
    % Arguments:
    %   ShuffleTestResults: one element of the struct array returned by
    %                        ShuffleTest_SpeechEEGCoherence.m (so ONE subject's results). Must
    %                        have fields True_MSC, Null_MSC, Chanlocs.
    %
    % Optional name-value arguments:
    %   options.Threshold (default norminv(0.95) ~= 1.6449): the per-channel rank-based-Z cutoff
    %       used to decide which channels are cluster CANDIDATES before spatial grouping. This is
    %       NOT itself a corrected significance threshold -- it only controls how generously
    %       channels are allowed to join a candidate cluster; the actual significance correction
    %       happens at the cluster-mass level (options.ClusterAlpha). Lower it to grow more
    %       permissive/larger candidate clusters, raise it to require each contributing channel
    %       to individually look stronger.
    %   options.ClusterAlpha (default 0.05): a real cluster is called significant if its
    %       cluster-mass p-value (see above) is <= this value.
    %   options.Adjacency (default []): pass a precomputed adjacency matrix (e.g. from
    %       BuildChannelAdjacency.m) to reuse across subjects/calls instead of rebuilding it from
    %       Chanlocs every time (they should all share the same montage, so this saves repeated
    %       work and guarantees every subject uses an IDENTICAL adjacency structure). If empty,
    %       it's built internally from ShuffleTestResults.Chanlocs with default settings.
    %
    % Returns ClusterResults, a struct with fields:
    %   .ChannelRankZ         [1 x Num_channels] the per-channel rank-based Z the whole test is
    %                         built on (identical to calling ComputeRankBasedZScore.m directly
    %                         on this subject's True_MSC/Null_MSC).
    %   .Clusters             {1 x Num_clusters} cell array; Clusters{k} = row vector of channel
    %                         INDICES making up real cluster k (channels with Z > Threshold that
    %                         are spatially connected).
    %   .ClusterLabels        {1 x Num_clusters} cell array of channel-label cell arrays (same
    %                         grouping as .Clusters, but as chanlocs(...).labels strings, for
    %                         readability).
    %   .ClusterMass          [1 x Num_clusters] sum of ChannelRankZ within each cluster.
    %   .ClusterSize          [1 x Num_clusters] number of channels in each cluster.
    %   .ClusterPeakZ         [1 x Num_clusters] max ChannelRankZ within each cluster.
    %   .ClusterPValue        [1 x Num_clusters] cluster-mass p-value (see above).
    %   .ClusterSignificant   [1 x Num_clusters] logical, ClusterPValue <= options.ClusterAlpha.
    %   .NullMaxClusterMass   [1 x NumShuffles] the max-cluster-mass null distribution itself
    %                         (useful for a sanity-check histogram against the real cluster
    %                         masses).
    %   .Threshold, .ClusterAlpha, .NumShuffles, .Adjacency: echoed-back settings.
    %   .Chanlocs, .Subj_id (if present on the input): passed through for convenience.
    %
    % Clusters are returned in no particular order; sort by .ClusterMass yourself if you want
    % "biggest cluster first" (e.g. [~, order] = sort(ClusterResults.ClusterMass, 'descend')).
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Why leave-one-out for a null draw's own Z, rather than scoring it against the full
    %    Null_MSC (including that draw) the way the real True_MSC's Z is scored: the real
    %    True_MSC is not one of the 1000 null draws, so no leave-one-out is needed for it. But a
    %    null draw IS one of the 1000 values in its own column of Null_MSC -- scoring it against
    %    a null distribution that includes itself would trivially and identically bias every
    %    null draw's Z upward relative to a hypothetical fresh draw, undermining the very null
    %    distribution this procedure is trying to build. Excluding it keeps every null draw's Z
    %    computed exactly the way the real data's Z is computed: against a null distribution it
    %    is not part of.
    %   2) This function operates on ONE subject at a time, by design -- see
    %    GroupLevel_ShuffleTest_Analysis.m for the (separate, still-random-effects) question of
    %    combining RankZ across subjects at FIXED channels. Combining cluster-level results
    %    ACROSS subjects, when the clusters may sit at different channels per subject, is a
    %    distinct group-level design question (e.g. population-prevalence-style inference) that
    %    is deliberately NOT decided inside this file.
    %   3) Not run in an actual MATLAB session (none was reachable while writing this) --
    %    syntax-checked only. Try it on one subject's ShuffleTestResults and sanity-check the
    %    output (cluster sizes, p-values, NullMaxClusterMass histogram) before trusting it.
    % =========================================================================================
    arguments
        ShuffleTestResults (1,1) struct
        options.Threshold (1,1) double {mustBeReal} = norminv(0.95)
        options.ClusterAlpha (1,1) double {mustBeReal, mustBePositive} = 0.05
        options.Adjacency (:,:) logical = logical([])
    end

    required_fields = {'True_MSC', 'Null_MSC', 'Chanlocs'};
    for f = 1:numel(required_fields)
        if ~isfield(ShuffleTestResults, required_fields{f})
            error('ClusterPermutationTest_ShuffleTest:MissingField', ...
                'ShuffleTestResults is missing required field ''%s''.', required_fields{f});
        end
    end

    True_MSC = reshape(ShuffleTestResults.True_MSC, 1, []);
    Null_MSC = ShuffleTestResults.Null_MSC;
    Chanlocs = ShuffleTestResults.Chanlocs;
    Num_channels = numel(True_MSC);
    NumShuffles  = size(Null_MSC, 1);

    if size(Null_MSC, 2) ~= Num_channels
        error('ClusterPermutationTest_ShuffleTest:SizeMismatch', ...
            'True_MSC has %d channels but Null_MSC has %d columns.', Num_channels, size(Null_MSC,2));
    end
    if numel(Chanlocs) ~= Num_channels
        error('ClusterPermutationTest_ShuffleTest:ChanlocsMismatch', ...
            'True_MSC has %d channels but Chanlocs has %d entries.', Num_channels, numel(Chanlocs));
    end

    % ---- Adjacency structure ------------------------------------------------------------------
    if isempty(options.Adjacency)
        [Adjacency, ~] = BuildChannelAdjacency(Chanlocs);
    else
        Adjacency = options.Adjacency;
        if ~isequal(size(Adjacency), [Num_channels, Num_channels])
            error('ClusterPermutationTest_ShuffleTest:AdjacencySizeMismatch', ...
                'options.Adjacency is %dx%d but there are %d channels.', ...
                size(Adjacency,1), size(Adjacency,2), Num_channels);
        end
    end

    % ---- Real data: rank-based Z, threshold, cluster -------------------------------------------
    ChannelRankZ = ComputeRankBasedZScore(True_MSC, Null_MSC);
    [real_clusters, real_mass, real_peak] = find_clusters(ChannelRankZ, options.Threshold, Adjacency);

    % ---- Null distribution of the max cluster mass, via leave-one-out on each shuffle draw -----
    NullMaxClusterMass = zeros(1, NumShuffles);
    all_idx = 1:NumShuffles;
    for s = 1:NumShuffles
        surrogate_true = Null_MSC(s, :);
        surrogate_null = Null_MSC(all_idx ~= s, :);
        surrogate_Z = ComputeRankBasedZScore(surrogate_true, surrogate_null);
        [~, s_mass, ~] = find_clusters(surrogate_Z, options.Threshold, Adjacency);
        if isempty(s_mass)
            NullMaxClusterMass(s) = 0;
        else
            NullMaxClusterMass(s) = max(s_mass);
        end
    end

    % ---- Cluster-level p-values (Phipson & Smyth, 2010 estimator, as used elsewhere in this
    %      codebase -- e.g. ShuffleTest_SpeechEEGCoherence.m's own p-value formula) -------------
    Num_clusters = numel(real_clusters);
    ClusterPValue = nan(1, Num_clusters);
    for k = 1:Num_clusters
        exceed_count = sum(NullMaxClusterMass >= real_mass(k));
        ClusterPValue(k) = (1 + exceed_count) / (NumShuffles + 1);
    end
    ClusterSignificant = ClusterPValue <= options.ClusterAlpha;

    ClusterLabels = cell(1, Num_clusters);
    for k = 1:Num_clusters
        ClusterLabels{k} = {Chanlocs(real_clusters{k}).labels};
    end

    ClusterResults = struct();
    ClusterResults.ChannelRankZ       = ChannelRankZ;
    ClusterResults.Clusters           = real_clusters;
    ClusterResults.ClusterLabels      = ClusterLabels;
    ClusterResults.ClusterMass        = real_mass;
    ClusterResults.ClusterSize        = cellfun(@numel, real_clusters);
    ClusterResults.ClusterPeakZ       = real_peak;
    ClusterResults.ClusterPValue      = ClusterPValue;
    ClusterResults.ClusterSignificant = ClusterSignificant;
    ClusterResults.NullMaxClusterMass = NullMaxClusterMass;
    ClusterResults.Threshold          = options.Threshold;
    ClusterResults.ClusterAlpha       = options.ClusterAlpha;
    ClusterResults.NumShuffles        = NumShuffles;
    ClusterResults.Adjacency          = Adjacency;
    ClusterResults.Chanlocs           = Chanlocs;
    if isfield(ShuffleTestResults, 'Subj_id')
        ClusterResults.Subj_id = ShuffleTestResults.Subj_id;
    end
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function [clusters, cluster_mass, cluster_peak] = find_clusters(Z, threshold, Adjacency)
    % Groups channels with Z > threshold into spatially contiguous clusters using Adjacency
    % (logical Num_channels x Num_channels), via a breadth-first search restricted to candidate
    % channels. Returns clusters (cell array of channel-index row vectors), cluster_mass (sum of
    % Z within each cluster), cluster_peak (max Z within each cluster). All three come back empty
    % (1x0 double / 1x0 cell) if no channel exceeds threshold.
    Num_channels = numel(Z);
    candidate = Z > threshold;
    visited = false(1, Num_channels);
    clusters = {};
    cluster_mass = [];
    cluster_peak = [];

    for ch = 1:Num_channels
        if candidate(ch) && ~visited(ch)
            % BFS over candidate, adjacent channels starting from ch.
            queue = ch;
            visited(ch) = true;
            this_cluster = ch;
            qi = 1;
            while qi <= numel(queue)
                current = queue(qi);
                qi = qi + 1;
                neighbors = find(Adjacency(current, :));
                for n = neighbors
                    if candidate(n) && ~visited(n)
                        visited(n) = true;
                        queue(end+1) = n; %#ok<AGROW>
                        this_cluster(end+1) = n; %#ok<AGROW>
                    end
                end
            end
            clusters{end+1} = sort(this_cluster); %#ok<AGROW>
            cluster_mass(end+1) = sum(Z(this_cluster)); %#ok<AGROW>
            cluster_peak(end+1) = max(Z(this_cluster)); %#ok<AGROW>
        end
    end
end
