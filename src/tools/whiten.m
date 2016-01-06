function X = whiten(X)
    X = X - min(X(:));
    X = X ./ max(X(:));