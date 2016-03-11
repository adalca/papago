function sm = model5exp(s, X, W, muk, dfn) 

    sm = 0;
    for i = 1:size(X, 1)
        x = X(i, :);
        w = W(i, :);
        d = dfn(w); 
        
        smtmp = -0.5 * logdet(s + d) - 0.5 * (x(:) - muk(:))' / (s + d) * (x(:) - muk(:));
        
        sm = sm + smtmp;
    end