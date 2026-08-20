function list_scrambled = derangeArray(original_list)
% DERANGEARRAY Randomly permute an array such that no element remains
% in its original position (a "derangement").
%
%   list_scrambled = derangeArray(original_list)

    N = numel(original_list);

    if N < 2
        error('Array must have at least 2 elements to create a derangement.');
    end

    original_list = original_list(:)';  % ensure row vector for indexing consistency
    origIdx = 1:N;

    % Rejection sampling: keep generating random permutations of the
    % indices until none of them map back to their original position.
    % For reasonably sized N this converges quickly (expected ~e ≈ 2.72
    % attempts), since P(no fixed points) -> 1/e as N grows.
    isDerangement = false;
    while ~isDerangement
        newIdx = randperm(N);
        isDerangement = all(newIdx ~= origIdx);
    end
    list_scrambled = original_list(newIdx);
end