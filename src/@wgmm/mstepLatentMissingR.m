function params = mstepLatentMissingR(wg, data)

    % data
    R = data.R;                     % R rotations from atlas to subject space.
    G = data.G;                     % Gamma rotations from subject space to atlas space
    Y = data.Y;                     % data in subject space
    ydsmasks = data.ydsmasks;       % down-sample mask cell -- (should be "planes" in subject space)
    
    % dimensions and numbers
    N = numel(ydsmasks);            % number of elements
    dHigh = size(R.data, 2);        % atlas-space dimension
    K = data.K;                     % number of clusters
    Nk = sum(wg.expect.gammank, 1); % sum of cluster support

    % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows -- so we
    % transpose the R data matrix and index in columns instead of rows.
    RdataTrans = R.data';
    
    % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than indexing
    % into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times per iteration, it
    % ends up mattering for our iterations.
    sigmaPrevCell = dimsplit(3, wg.params.sigma);
    
    % update Yhat in atlas space. 
    %   This should really be part of the E-step, but we'll do it here to
    %   not pass large matrices around. Maybe we should fix this for clenliness!
    Yhat = zeros(N, dHigh, K);
    for i = 1:N
        % extract the observed y values in original space.
        obsIdx = ydsmasks{i};
        ySubj = Y{i};
        ySubjObs = ySubj(obsIdx);
        
        % prepare the rotation matrices
        r = RdataTrans(:, R.idx{i})';
        robs = r(obsIdx, :);
        rmis = r(~obsIdx, :);
        g = G.data(:, G.idx{i}); % pinv(r) might actually be better. see testGammaRotationMatrices().

        % go through each cluster
        for k = 1:K
            % extract the kth statistics
            sigma = sigmaPrevCell{k};
            mu = wg.params.mu(k, :);
            
            % rotate statistics to subject space. see paffine.atl2SubjGauss().
            %   doing sigmaSubj = r * (sigma * robs');: since the right multiplication gives fewer entries, this is much faster.
            % update estimate missing values
            srobs = (sigma * robs');
            oosigmak = robs * srobs;
            mosigmak = rmis * srobs;
            muSubj = mu * r';

            ySubjk = ySubj(:)';
            ySubjk(~obsIdx) = muSubj(~obsIdx) + (mosigmak * (oosigmak \ (ySubjObs - muSubj(obsIdx))'))';

            % rotate back to atlas space
            Yhat(i, :, k) = ySubjk * g';
        end
    end
    
    % mu update (same as in latentMissing, using Yhat)
    params.mu = zeros(K, dHigh);
    for k = 1:K
        gnk = wg.expect.gammank(:, k);
        params.mu(k, :) = sum(bsxfun(@times, gnk, Yhat(:, :, k))) ./ Nk(k);
    end
    
    % update new sigma using new mu
    sigma = arrayfunc(@(s) zeros(dHigh, dHigh), 1:K);
    for i = 1:N
        % extract the observed y values in original space.
        obsIdx = ydsmasks{i};
        
        % prepare the rotation matrices
        r = RdataTrans(:, R.idx{i})';
        robs = r(obsIdx, :);
        rmis = r(~obsIdx, :);
        
        % gamma rotation for missing values
        %   gMissing rows won't add up to 1, but we only multiply gMissing with a sparse S which
        %   would have 0s in the gObserved entries, so the result should be the same
        gMissing = G.data(:, G.idx{i}(~obsIdx));
        
        % go through each cluster
        for k = 1:K
            gnk = wg.expect.gammank(i, k);
            
            % compute S correction term
            sigmaPrev = sigmaPrevCell{k};
            
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
             
            % Method 2: only compute what we need in Subject space. This saves a bit of time (5%?)
            % extract appropriate sigmas
            sr = (sigmaPrev * robs');
            oosigmak = robs * sr;
            mosigmak = rmis * sr;
            % compute correction term
            sCorrTerm1 = rmis * (sigmaPrev * rmis');
            sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
            sCorrSubjSpace = sCorrTerm1 - sCorrTerm2;
            sCorrAtlSpace = gMissing * sCorrSubjSpace * gMissing';       

            % empirical covariance for this data point
            ySubj = Yhat(i, :, k);
            sEmpirical = (ySubj - params.mu(k,:))' * (ySubj - params.mu(k,:)); 

            % update sigma
            sigma{k} = sigma{k} + gnk *  (sCorrAtlSpace + sEmpirical);
        end
    end
    % combine sigma to [nData-by-nData-by-K]
    params.sigma = cat(3, sigma{:});
    assert(isclean(params.sigma), 'sigma is not clean');

    % normalize sigma
    for k = 1:K
        params.sigma(:, :, k) = params.sigma(:, :, k) ./ Nk(k);
    end
    
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
