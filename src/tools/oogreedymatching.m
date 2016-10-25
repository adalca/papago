function c = oogreedymatching(V1, V2)
% vectors are rows
% so V1 is Nvecs * D

    c = zeros(1, size(V1, 1));

    for i = 1:size(V1, 1)
        [~, mi] = min(pdist2(V1(i, :), V2));
        c(i) = mi;
        V2(mi, :) = inf;
    end
    