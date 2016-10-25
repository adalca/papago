function r = sigmaUpdateEqn(X, s, mu, W, w2Dfn)

    r = 0;
    for i = 1:size(X, 1)
        x = X(i, :);
        w = W(i, :);
        D = w2Dfn(w);
        
        isd = inv(s + D);
        xmd = (x - mu);
        dfd = xmd(:) * xmd(:)';
        
        r = s + isd - isd * dfd * isd;
    end
    