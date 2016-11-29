function params = mstepLatentMissingR(wg, data)

    % data
    R = data.R;                     % R rotations from atlas to subject space.
    G = data.G;                     % Gamma rotations from subject space to atlas space
    Y = data.Y;                     % data in subject space
    ydsmasks = data.ydsmasks;       % down-sample mask cell -- (should be "planes" in subject space)
    ydsmasksFullVoxels = data.ydsmasksFullVoxels;
    
    % dimensions and numbers
    N = numel(ydsmasks);            % number of elements
    dHigh = size(R.data, 2);        % atlas-space dimension
    K = data.K;                     % number of clusters
    Nk = sum(wg.expect.gammank, 1); % sum of cluster support

    % throw a warning reminding us of a bunch of choices we're making about how to deal with data on
    % the edges. First, we only re-estimate data in the subject space that we can contribute to with
    % more than misThr weight of the atlas space voxels. Then we only use those voxels (and the
    % observed voxels) in computing Yhat (the re-estimated Y in atlas space).
    warning('Choice of which voxels to include in rotation are uncertain.');
    misThr = 0.95;
    % obsThr = 0.97; % should use for obsIdx?    
    
    % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows -- so we
    % transpose the R data matrix and index in columns instead of rows.
    RdataTrans = R.data';
    
    % prepare new variables
    newMu = zeros(K, dHigh);
    newSigma = zeros(dHigh, dHigh, K);
    
    % Operate per cluster since there is no interaction among clusters, and this way we can keep
    % much smaller variables around.
    for k = 1:K
        sigmaPrev = wg.params.sigma(:, :, k);
        muPrev = wg.params.mu(k, :);
        expectk = wg.expect.gammank(:, k);

        % update Yhat in atlas space. 
        %   This should really be part of the E-step, but we'll do it here to
        %   not pass large matrices around.
        Yhat = zeros(N, dHigh);
        for i = 1:N
            % don't need to compute Yhat_ki, since it's only used when gamma_ki is > 0
            % TODO: should make sure to include this in the normalization.
            if expectk(i) < 1e-4, continue; end
            
            % extract the observed y values in original space.
            obsIdx = ydsmasksFullVoxels{i};
            misIdx = data.rWeight{i} > misThr & ~ydsmasks{i}; % used to be just ~obsIdx
            ySubj = Y{i};
            ySubjObs = ySubj(obsIdx);

            % prepare the rotation matrices
            r = RdataTrans(:, R.idx{i})';
            robs = r(obsIdx, :);
            rmis = r(misIdx, :);

            % rotate statistics to subject space. see paffine.atl2SubjGauss().
            %   It's important to pay attention to order: 
            %   >> sigmaSubj = r * (sigma * robs'); 
            %   is fast since the right multiplication gives fewer entries
            srobs = (sigmaPrev * robs');
            oosigmak = robs * srobs;
            mosigmak = rmis * srobs;
            muSubj = muPrev * r';

            % update Y
            ySubjk = ySubj(:)';
            ySubjk(misIdx) = muSubj(misIdx) + (mosigmak * (oosigmak \ (ySubjObs - muSubj(obsIdx))'))';

            % rotate back to atlas space
            %   - pinv(r) might actually be better. see testGammaRotationMatrices()
            %   - need to avoid values we didn't re-estimate and were not observed
            validVals = misIdx | ydsmasks{i};
            g = G.data(:, G.idx{i}(validVals)); % g is atl-by-subj
            g = bsxfun(@rdivide, g, sum(g, 2)); % make sure it sums to 1.
            Yhat(i, :) = ySubjk(validVals) * g';
        end
        f = isnan(Yhat);
        if sum(f(:)) > 0
            warning('mstepLMR: Found %d NANs in Yhat. Filling them in with the mean', sum(f(:)));
            [~, b] = ind2sub(size(Yhat), find(f));
            Yhat(f) = muPrev(b);
        end
        assert(isclean(Yhat));

        % mu update (same as in latentMissing, using Yhat)
        newMuk = sum(bsxfun(@times, expectk, Yhat)) ./ sum(expectk);

        % update sigma. Start with the empirical covariance. Note:
        %   - ECMNMLE use the most recent mu to *recompute* y_ik. This does not make sense to us
        %     y_ik = Y(i, :);
        %     y_ik(~obsIdx) = newMuk(~obsIdx) + (mosigmak * (oosigmak \ (Y(i, obsIdx) - newMuk(obsIdx))'))';
        ydiff = bsxfun(@times, sqrt(wg.expect.gammank(:, k)), bsxfun(@minus, Yhat, newMuk));
        newSigmak = ydiff' * ydiff;
        
        for i = 1:N
            if expectk(i) < 1e-4, continue; end
            
            % extract the observed y values in original space.
            obsIdx = ydsmasksFullVoxels{i};
            misIdx = data.rWeight{i} > misThr & ~ydsmasks{i}; % used to be ~obsIdx        

            % prepare the rotation matrices
            r = RdataTrans(:, R.idx{i})';
            robs = r(obsIdx, :);
            rmis = r(misIdx, :);

            % gamma rotation for missing values
            %   gMissing rows won't add up to 1, but we only multiply gMissing with a sparse S which
            %   would have 0s in the gObserved entries, so the result should be the same
            gMissing = G.data(:, G.idx{i}(misIdx));

            % only compute what we need in Subject space. This saves a bit of time (5%?)
            % extract appropriate sigmas
            sr = (sigmaPrev * robs');
            oosigmak = robs * sr;
            mosigmak = rmis * sr;
            % compute correction term
            sCorrTerm1 = rmis * (sigmaPrev * rmis');
            sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
            sCorrSubjSpace = sCorrTerm1 - sCorrTerm2;
            sCorrAtlSpace = gMissing * sCorrSubjSpace * gMissing';       

            % update sigma
            newSigmak = newSigmak + expectk(i) *  sCorrAtlSpace;
        end
        
        % normalize and record stats
        newMu(k, :) = newMuk;
        newSigma(:, :, k) = newSigmak ./ sum(expectk);
    end
        
    % update parameters
    params.mu = newMu;
    params.sigma = newSigma;
    assert(isclean(params.mu), 'mu is not clean');
    assert(isclean(params.sigma), 'sigma is not clean');
    
    % estimate low dimensional representation
    %   estimate sigmasq, W (principal components)
    %   update sigma to be C = W'W + sI. 
    if wg.opts.model.dopca % need to also return W, sigmasq
        
        % compute the covariance of sigma.
        for k = 1:K
            % SVD
            [u, srobs, ~] = svd(params.sigma(:,:,k));
            ds = diag(srobs);
            
            % PCA components
            dLow = wg.opts.model.dopca;
            sigmasq = 1 ./ (dHigh - dLow) * sum(ds(dLow+1:end));
            W = u(:, 1:dLow) * sqrt(srobs(1:dLow, 1:dLow) - sigmasq * eye(dLow));
            Covar = W*W' + sigmasq*eye(size(u,2));

            % save results
            params.sigmasq(k) = sigmasq;
            params.W(:,:,k) = W;
            params.sigma(:,:,k) = Covar;
        end
    end
    
    % add regularization diagonal, check for PDness, etc.
    [params.sigma, params.sigmainv] = wgmm.sigmafull(wg, params.sigma);
    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));
end

% old sigmak computations
% Method 1: compute entire Subject-space Sigma, then take parts of it.
% % rotate sigma. see paffine.atl2SubjGauss
% sigmaPrevSubj = r * sigmaPrev * r'; 
% % extract appropriate sigmas
% oosigmak = sigmaPrevSubj(obsIdx, obsIdx);
% mosigmak = sigmaPrevSubj(~obsIdx, obsIdx);    
% % compute correction term
% sCorrTerm1 = sigmaPrevSubj(~obsIdx, ~obsIdx);
% sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
% sCorrSubjSpace = sCorrTerm1 - sCorrTerm2;
% sCorrAtlSpace = gMissing * sCorrSubjSpace * gMissing'; 
