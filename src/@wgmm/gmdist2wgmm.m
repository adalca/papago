function gmm = gmdist2wgmm(gmdist, X, W, K)
    gmsigmainv = zeros(size(gmdist.Sigma));
    for i = 1:size(gmdist.Sigma, 3)
        gmsigmainv(:,:,i) = inv(gmdist.Sigma(:,:,i));
    end

    gmm = wgmm(gmdist.mu, gmdist.Sigma, gmdist.ComponentProportion, gmsigmainv);
end