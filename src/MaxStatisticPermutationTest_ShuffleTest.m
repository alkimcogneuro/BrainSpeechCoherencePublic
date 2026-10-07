function [MaxStatResults] = MaxStatisticPermutationTest_ShuffleTest(ShuffleTestResults, options)
    % =========================================================================================
    % Single-Subject Max-Statistic Permutation Test for Shuffle-Test MSC Results (No Clustering)
    % =========================================================================================
    % A simpler, spatially-blind alternative to ClusterPermutationTest_ShuffleTest.m, added as an
    % independent sanity check while that file's channel-adjacency/clustering logic is still
    % being verified. This asks a narrower question with NO reference to electrode positions or
    % neighbors at all: is there at least one channel whose true MSC exceeds what the single BEST
    % channel achieves purely by chance?
    %
    % This is the same "family-wise error correction via the max statistic" idea used inside
    % ClusterPermutationTest_ShuffleTest.m's own null (see Nichols & Holmes, 2002, for this
    % approach in neuroimaging generally) -- it's exactly what that file's cluster test reduces
    % to if you fed it an adjacency matrix with NO edges at all (every channel its own cluster of
    % size 1), just computed here directly, without ever building or touching an adjacency
    % structure. If this test and the cluster test disagree in some obviously wrong way (e.g.
    % this finds a clearly strong single channel that the cluster test somehow misses or badly
    % mislocates), that points at the clustering step (BuildChannelAdjacency.m / the
    % connected-components logic) rather than at the underlying shuffle-test statistics, which
    % this function and the cluster test both start from identically.
    %
    % Procedure:
    %   1) options.Metric='MSC' (default): work directly in True_MSC / Null_MSC units, the same
    %      units ShuffleTest_SpeechEEGCoherence.m's own per-channel p-values already use -- no
    %      extra transform, so this is about as close to "just re-deriving the base shuffle
    %      test's logic, but corrected for scanning every channel" as possible, and easiest to
    %      hand-check against ShuffleTestResults.PValues/True_MSC directly.
    %      options.Metric='RankZ': operate on ComputeRankBasedZScore.m's rank-based Z instead
    %      (computed per-channel for the true data, and via leave-one-out for each null draw --
    %      same convention ClusterPermutationTest_ShuffleTest.m uses). More statistically
    %      appropriate if different channels' null distributions have noticeably different
    %      spread (raw-MSC max-statistic testing can be dominated by whichever channel happens to
    %      have the noisiest/highest-variance null, not necessarily the most genuinely coherent
    %      one) -- but 'MSC' is the more direct, easier-to-hand-check default for this pass.
    %   2) For the real data, take the per-channel statistic's max across channels
    %      (True_MaxStat / BestChannelStat below), and note which channel achieved it.
    %   3) For each of the NumShuffles null draws, take ITS max across channels too (using the
    %      same per-draw leave-one-out RankZ if options.Metric='RankZ'; raw MSC values need no
    %      leave-one-out since nothing is being ranked against anything for that metric). This
    %      gives NullMaxStat, a [1 x NumShuffles] null distribution of "the best channel's value
    %      under pure chance."
    %   4) BestChannelStat's p-value = (1 + count(NullMaxStat >= BestChannelStat)) /
    %      (NumShuffles + 1) (same Phipson & Smyth, 2010 estimator used everywhere else in this
    %      codebase). This single p-value already answers "is there at least one channel
    %      significantly better than chance, correcting for having looked at all of them."
    %   5) As a by-product of the same null, every INDIVIDUAL channel also gets its own
    %      max-corrected p-value ((1 + count(NullMaxStat >= that channel's stat)) / (NumShuffles
    %      + 1)) -- a full family-wise-error-corrected per-channel significance map, at no extra
    %      cost, and a legitimate (if conservative) alternative to per-channel FDR.
    %
    % Arguments:
    %   ShuffleTestResults: one element of the struct array returned by
    %                        ShuffleTest_SpeechEEGCoherence.m (one subject). Must have fields
    %                        True_MSC, Null_MSC (Chanlocs used only for channel labels, if
    %                        present).
    %
    % Optional name-value arguments:
    %   options.Metric (default 'MSC'): 'MSC' or 'RankZ' -- see above.
    %   options.Alpha   (default 0.05): significance threshold applied to both the overall
    %                   best-channel p-value and each channel's individual max-corrected p-value.
    %
    % Returns MaxStatResults, a struct with fields:
    %   .Metric                 which metric was used ('MSC' or 'RankZ').
    %   .ChannelStat             [1 x Num_channels] the per-channel statistic actually used.
    %   .ChannelPValue           [1 x Num_channels] each channel's max-corrected p-value.
    %   .ChannelSignificant      [1 x Num_channels] logical, ChannelPValue <= options.Alpha.
    %   .BestChannelIndex        index of the channel with the largest ChannelStat.
    %   .BestChannelLabel        that channel's label (from Chanlocs, if present).
    %   .BestChannelStat         ChannelStat at BestChannelIndex (== max(ChannelStat)).
    %   .BestChannelPValue       ChannelPValue at BestChannelIndex (== the overall test's p-value).
    %   .HasSignificantChannel   logical, true if BestChannelPValue <= options.Alpha (equivalently,
    %                            any(ChannelSignificant) -- same test, since the max-corrected
    %                            p-value is monotonic in ChannelStat).
    %   .NullMaxStat             [1 x NumShuffles] the max-over-channels null distribution itself.
    %   .Alpha, .NumShuffles: echoed-back settings.
    %   .Chanlocs, .Subj_id (if present on the input): passed through for convenience.
    %
    % Not run in an actual MATLAB session (none was reachable while writing this) --
    % syntax-checked only, same caveat as the rest of this analysis family.
    % =========================================================================================
    arguments
        ShuffleTestResults (1,1) struct
        options.Metric (1,:) char {mustBeMember(options.Metric, {'MSC', 'RankZ'})} = 'MSC'
        options.Alpha (1,1) double {mustBeReal, mustBePositive} = 0.05
    end

    required_fields = {'True_MSC', 'Null_MSC'};
    for f = 1:numel(required_fields)
        if ~isfield(ShuffleTestResults, required_fields{f})
            error('MaxStatisticPermutationTest_ShuffleTest:MissingField', ...
                'ShuffleTestResults is missing required field ''%s''.', required_fields{f});
        end
    end

    True_MSC = reshape(ShuffleTestResults.True_MSC, 1, []);
    Null_MSC = ShuffleTestResults.Null_MSC;
    Num_channels = numel(True_MSC);
    NumShuffles  = size(Null_MSC, 1);

    if size(Null_MSC, 2) ~= Num_channels
        error('MaxStatisticPermutationTest_ShuffleTest:SizeMismatch', ...
            'True_MSC has %d channels but Null_MSC has %d columns.', Num_channels, size(Null_MSC,2));
    end

    % ---- Per-channel statistic, for the real data and for every null draw ---------------------
    switch options.Metric
        case 'MSC'
            ChannelStat = True_MSC;               % [1 x Num_channels]
            NullStat = Null_MSC;                  % [NumShuffles x Num_channels]
        case 'RankZ'
            ChannelStat = ComputeRankBasedZScore(True_MSC, Null_MSC);
            NullStat = zeros(NumShuffles, Num_channels);
            all_idx = 1:NumShuffles;
            for s = 1:NumShuffles
                NullStat(s, :) = ComputeRankBasedZScore(Null_MSC(s, :), Null_MSC(all_idx ~= s, :));
            end
    end

    % ---- Null distribution of the max-across-channels statistic --------------------------------
    NullMaxStat = max(NullStat, [], 2);  % [NumShuffles x 1]

    % ---- Per-channel and overall (best-channel) max-corrected p-values -------------------------
    ChannelPValue = nan(1, Num_channels);
    for ch = 1:Num_channels
        ChannelPValue(ch) = (1 + sum(NullMaxStat >= ChannelStat(ch))) / (NumShuffles + 1);
    end
    ChannelSignificant = ChannelPValue <= options.Alpha;

    [BestChannelStat, BestChannelIndex] = max(ChannelStat);
    BestChannelPValue = ChannelPValue(BestChannelIndex);

    Chanlocs = [];
    if isfield(ShuffleTestResults, 'Chanlocs') && numel(ShuffleTestResults.Chanlocs) == Num_channels
        Chanlocs = ShuffleTestResults.Chanlocs;
        BestChannelLabel = Chanlocs(BestChannelIndex).labels;
    else
        BestChannelLabel = sprintf('Channel_%d', BestChannelIndex);
    end

    MaxStatResults = struct();
    MaxStatResults.Metric                = options.Metric;
    MaxStatResults.ChannelStat           = ChannelStat;
    MaxStatResults.ChannelPValue         = ChannelPValue;
    MaxStatResults.ChannelSignificant    = ChannelSignificant;
    MaxStatResults.BestChannelIndex      = BestChannelIndex;
    MaxStatResults.BestChannelLabel      = BestChannelLabel;
    MaxStatResults.BestChannelStat       = BestChannelStat;
    MaxStatResults.BestChannelPValue     = BestChannelPValue;
    MaxStatResults.HasSignificantChannel = BestChannelPValue <= options.Alpha;
    MaxStatResults.NullMaxStat           = reshape(NullMaxStat, 1, []);
    MaxStatResults.Alpha                 = options.Alpha;
    MaxStatResults.NumShuffles           = NumShuffles;
    if ~isempty(Chanlocs)
        MaxStatResults.Chanlocs = Chanlocs;
    end
    if isfield(ShuffleTestResults, 'Subj_id')
        MaxStatResults.Subj_id = ShuffleTestResults.Subj_id;
    end
end
