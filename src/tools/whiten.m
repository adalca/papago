function X = whiten(X)
% force 0 mean and max 1
    X = X - min(X(:));
    X = X ./ max(X(:));