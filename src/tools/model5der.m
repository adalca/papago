function sm = model5der(s, X, W, muk, gammank, dfn) 

    sm = 0;
    for i = 1:size(X, 1)
        x = X(i, :);
        xmm = x(:) - muk(:);
        w = W(i, :);
        d = dfn(w); 
        
        
        sd = s + d;
        isd = inv(sd);
        
        smtmp = isd - isd * xmm * xmm * isd; % could speed up
        
        sm = sm + gammank(i) * smtmp;
    end