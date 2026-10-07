function [RankZ, Percentile, ExceedCount] = ComputeRankBasedZScore(True_value, Null_values)
    % =========================================================================================
    % Rank-Based (Percentile) Z-Score: A Distribution-Free Standardization of an Observed Value
    % Against an Empirical Null Distribution
    % =========================================================================================
    % Converts an observed value (e.g. ShuffleTestResults.True_MSC) into a signed z-like score
    % using its PERCENTILE RANK within an empirical null distribution (e.g.
    % ShuffleTestResults.Null_MSC), rather than the null's raw mean and standard deviation.
    %
    % Why this exists, instead of just using (True_value - mean(Null_values)) / std(Null_values)
    % (which is what ShuffleTestResults.ZScores computes): that raw standardization implicitly
    % treats the null distribution as roughly symmetric/Gaussian-shaped, which band-averaged MSC
    % is not guaranteed to be -- MSC is non-negative by construction (bounded at 0), and coherence
    % estimates are classically right-skewed, especially with modest trial counts. That skew and
    % hard floor mean the raw Z-score's negative tail is mechanically compressed relative to its
    % positive tail (True_MSC can only fall as low as 0, no matter how large Null_Std is), so equal
    % positive and negative raw Z-scores are NOT equally-extreme evidence, even though a symmetric
    % colormap or a naive reading of "Z" implies they are.
    %
    % This function sidesteps that entirely by never looking at the null distribution's shape,
    % mean, or standard deviation -- only at WHERE True_value RANKS among the null draws. Percentile
    % rank is, by construction, uniform on (0,1) under the null regardless of the underlying
    % distribution's shape (that's the core premise of a permutation test). Passing that percentile
    % through the inverse standard normal CDF (norminv) then gives a value that behaves like a
    % standard normal deviate under the null -- REGARDLESS of whether the raw values (True_MSC,
    % Null_MSC) are themselves Gaussian, skewed, or bounded. This is the same "probit transform of
    % a p-value" trick used in Stouffer's method for combining independent p-values.
    %
    % A continuity correction (rank - 0.5, rather than rank) is used so the output is always
    % finite, even when True_value beats or loses to every single null draw -- with the plain
    % (1+exceed_count)/(N+1) p-value formula, the "loses to everyone" case gives p exactly 1, which
    % maps to norminv(0) = -Inf without this correction.
    %
    % Arguments:
    %   True_value:  [1 x Num_channels] (or [Num_channels x 1]) observed value, one per channel
    %                (e.g. ShuffleTestResults.True_MSC).
    %   Null_values: [NumShuffles x Num_channels] empirical null draws, one column per channel,
    %                matching True_value's channel order (e.g. ShuffleTestResults.Null_MSC).
    %
    % Returns:
    %   RankZ:        [1 x Num_channels] the signed rank-based z-score. Large positive = True_value
    %                 sits near the top of its null distribution (strong evidence of an
    %                 above-chance effect). Large negative = True_value sits near the bottom of its
    %                 null distribution (True_value was unusually LOW relative to chance -- worth
    %                 treating as a data-quality flag for that channel/subject, not as evidence of
    %                 a below-chance effect in the opposite-but-equal sense a symmetric number
    %                 might suggest; see ShuffleTest_SpeechEEGCoherence.m's one-sided framing).
    %                 Near 0 = True_value looks like a typical/unremarkable draw from the null.
    %   Percentile:   [1 x Num_channels] the continuity-corrected percentile rank used to compute
    %                 RankZ (in (0,1), never exactly 0 or 1). RankZ = norminv(Percentile).
    %   ExceedCount:  [1 x Num_channels] number of null draws >= True_value, per channel -- the
    %                 same quantity ShuffleTest_SpeechEEGCoherence.m uses for its own p-value, kept
    %                 here so you can cross-check PValues = (1+ExceedCount)/(NumShuffles+1) against
    %                 whatever ShuffleTestResults.PValues already has for the same data.
    %
    % -----------------------------------------------------------------------------------------
    % Resolution caveat
    % -----------------------------------------------------------------------------------------
    % Like the underlying permutation p-value, RankZ's resolution is limited by NumShuffles: with
    % 1000 shuffles, "True_value beat all 1000 draws" and "True_value would have beaten 999,999 out
    % of a million draws" are indistinguishable -- both map to the same (near-)extreme RankZ. Very
    % large |RankZ| should be read as "at the resolution limit of this shuffle count," not as an
    % arbitrarily precise magnitude.
    %
    % No shuffles need to be re-run to use this: it only needs the True_value and Null_values
    % ShuffleTest_SpeechEEGCoherence.m already computed and stored.
    % =========================================================================================
    arguments
        True_value (1,:) double
        Null_values (:,:) double
    end

    True_value = reshape(True_value, 1, []);  % [1 x Num_channels]
    Num_channels = numel(True_value);

    if size(Null_values, 2) ~= Num_channels
        error('ComputeRankBasedZScore:ChannelMismatch', ...
            ['True_value has %d channels but Null_values has %d columns. Null_values must be ' ...
             '[NumShuffles x Num_channels], with columns matching True_value''s channel order ' ...
             '(e.g. ShuffleTestResults.True_MSC and ShuffleTestResults.Null_MSC).'], ...
            Num_channels, size(Null_values, 2));
    end

    NumShuffles = size(Null_values, 1);
    if NumShuffles < 1
        error('ComputeRankBasedZScore:NoNullDraws', 'Null_values has zero rows -- no null draws to rank against.');
    end

    % ---- Rank True_value within its own null distribution, per channel -----------------------
    % Same ">=" convention ShuffleTest_SpeechEEGCoherence.m uses for its own p-value, so
    % ExceedCount here is directly comparable to (and can be used to recompute/cross-check) that
    % file's PValues.
    ExceedCount = sum(Null_values >= True_value, 1);  % [1 x Num_channels], True_value broadcasts across rows

    % ---- Continuity-corrected percentile and its probit (inverse normal CDF) transform --------
    % p_cc is the one-sided permutation p-value with a continuity correction (0.5 instead of 1 in
    % the numerator) so it never hits exactly 0 or 1 for any finite NumShuffles, which keeps
    % norminv finite at both ends. p_cc ranges over roughly (0, 1), symmetric about 0.5 when
    % ExceedCount is at the midpoint of its possible range.
    p_cc = (0.5 + ExceedCount) / (NumShuffles + 1);
    Percentile = 1 - p_cc;
    RankZ = norminv(Percentile);

    % ---- Flag (but do not hide) non-finite results ---------------------------------------------
    % Should not occur given the continuity correction above (p_cc is always strictly inside
    % (0,1) for finite NumShuffles), but guarded explicitly in case Null_values contains NaN/Inf
    % entries that corrupted the ExceedCount computation upstream.
    bad_channels = find(~isfinite(RankZ));
    if ~isempty(bad_channels)
        warning('ComputeRankBasedZScore:NonFiniteResult', ...
            ['RankZ is non-finite for %d channel(s): %s. This likely means Null_values or ' ...
             'True_value contained NaN/Inf for those channels -- check the upstream ' ...
             'ShuffleTestResults for that subject.'], numel(bad_channels), mat2str(bad_channels(:)'));
    end
end
