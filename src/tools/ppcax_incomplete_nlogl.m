function nloglk_new = ppcax_incomplete_nlogl(Y, mu, Wnew, vnew)
% Y is p-by-n
% mu is a vector of size p

    obs = ~isnan(Y);
    mu = mu(:);

%     Compute negative log-likelihood function
    nloglk_new = 0;
    for j = 1:size(obs, 2)
        idxObs = obs(:,j);
        y = Y(idxObs,j) - mu(idxObs); % the jth observation centered with only complete elements
        Wobs = Wnew(idxObs,:);
        Cy = Wobs*Wobs'+vnew*eye(sum(idxObs));

        tt1 = sum(idxObs)*log(2*pi);
        tt2 = logdet(Cy); % AVD: logdet function instead of log(det(.))
        tt3 = trace(Cy\(y*y'));
        nloglk_new = nloglk_new + ...
            (tt1 + tt2 + tt3)/2;
    end
    assert(isclean(nloglk_new)); % avd addition