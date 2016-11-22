function params = mstepLatentMissingR(wg, data)

    % data
    R = data.R;                     % R rotations from atlas to subject space.
    G = data.G;                     % Gamma rotations from subject space to atlas space
    Y = data.Y;                     % data in subject space
    ydsmasks = data.ydsmasks;       % down-sample mask cell -- (should be "planes" in subject space)
    
    % dimensions and numbers
    N = numel(ydsmasks);            % number of elements
    dHigh = size(R.data, 2);        % atlas-space dimension
    K = size(wg.expect.gammank, 2); % number of clusters
    Nk = sum(wg.expect.gammank, 1); % sum of cluster support

    % Prepare existing sigma as K cells. Indexing into a small cell is faster than indexing into a
    % 3-D matrix [nData-by-nData-by-K], and since we access it this K times per iteration, it ends
    % up mattering for our iterations.
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
        r = R.data(R.idx{i}, :);
        g = G.data(:, G.idx{i}); % pinv(r) might actually be better. see testGammaRotationMatrices().

        % go through each cluster.
        for k = 1:K
            % extract the kth statistics
            sigma = sigmaPrevCell{k};
            mu = wg.params.mu(k, :);

            % rotate statistics to subject space. see paffine.atl2SubjGauss().
            sigmaSubj = r * sigma * r(obsIdx, :)';
            muSubj = mu * r';

            % update estimate missing values
            oosigmak = sigmaSubj(obsIdx, :);
            mosigmak = sigmaSubj(~obsIdx, :);

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
        params.mu(k, :) = sum(bsxfun(@times, gnk, Yhat(:, :, k))) ./ sum(gnk);
    end
    
    % update new sigma using new mu
    sigma = arrayfunc(@(s) zeros(dHigh, dHigh), 1:K);
    for i = 1:N
        % extract the observed y values in original space.
        obsIdx = ydsmasks{i};
        
        % prepare the rotation matrices
        r = R.data(R.idx{i}, :);
        g = G.data(:, G.idx{i});

        for k = 1:K
            gnk = wg.expect.gammank(i, k);
            
            % gamma rotation for missing values
            %   gMissing rows won't add up to 1, but we only multiply gMissing with a sparse S which
            %   would have 0s in the gObserved entries, so the result should be the same
            gMissing = g(:, ~obsIdx); 
            
            % compute S correction term
            sigmaPrev = sigmaPrevCell{k};
            % rotate sigma. see paffine.atl2SubjGauss
            sigmaPrevSubj = r * sigmaPrev * r';
            % extract appropriate sigmas
            oosigmak = sigmaPrevSubj(obsIdx, obsIdx);
            mosigmak = sigmaPrevSubj(~obsIdx, obsIdx);    
            % compute correction term
            sCorrTerm1 = sigmaPrevSubj(~obsIdx, ~obsIdx);
            sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
            sCorrSubjSpace = sCorrTerm1 - sCorrTerm2;
            sCorrAtlSpace = gMissing * sCorrSubjSpace * gMissing'; 

            % empirical covariance for this data point
            y_ik = Yhat(i, :, k);
            sEmpirical = (y_ik - params.mu(k,:))' * (y_ik - params.mu(k,:)); 

            % update sigma
            sigma{k} = sigma{k} + gnk *  (sCorrAtlSpace + sEmpirical);
        end
    end
    params.sigma = cat(3, sigma{:});
    assert(isclean(params.sigma), 'sigma is not clean');

    % estimate sigmasq, W (principal components)
    % update sigma to be C = W'W + sI. 
    if wg.opts.model.dopca % need to also return W, sigmasq
        
        % compute the covariance of sigma.
        for k = 1:K
            % SVD
            [u, s, ~] = svd(params.sigma(:,:,k));
            ds = diag(s);
            
            % PCA components
            dLow = wg.opts.model.dopca;
            sigmasq = 1 ./ (dHigh - dLow) * sum(ds(dLow+1:end));
            wts = u(:, 1:dLow) * sqrt(s(1:dLow, 1:dLow) - sigmasq * eye(dLow));
            Covar = wts*wts' + sigmasq*eye(size(u,2));

            params.sigmasq(k) = sigmasq;
            params.W(:,:,k) = wts;
            params.sigma(:,:,k) = Covar;
        end
    end
    
    % normalize
    for k = 1:K
        params.sigma(:, :, k) = params.sigma(:, :, k) ./ sum(wg.expect.gammank(:, k));
    end
    
    [params.sigma, params.sigmainv] = wgmm.sigmafull(wg, params.sigma);
    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));
end
