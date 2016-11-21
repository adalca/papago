function params = mstepLatentMissingR(wg, data)

    % data
    R = data.R;
    G = data.G;
    ydsmasks = data.ydsmasks;
    Y = data.Y;
    N = numel(ydsmasks);
    dHigh = size(R.data, 2);
    K = size(wg.expect.gammank, 2);
    Nk = sum(wg.expect.gammank, 1);

    % update Yhat in atlas space. This should really be part of the E-step, but we'll do it here to
    % not pass large matrices around. Maybe we should fix this!
    Yhat = zeros(N, dHigh, K);
    for i = 1:N
        obsIdx = ydsmasks{i};
        yobs = Y{i}(obsIdx);
        r = R.data(R.idx{i}, :);
        g = G.data(:, G.idx{i});

        for k = 1:K
            sigma = wg.params.sigma(:, :, k);
            mu = wg.params.mu(k, :);

            % rotate sigma. see paffine.atl2SubjGauss
            sigmaSubj = r * sigma * r(obsIdx, :)';
            muSubj = mu * r';

            % update Yhat_ik
            oosigmak = sigmaSubj(obsIdx, :);
            mosigmak = sigmaSubj(~obsIdx, :);

            y_ijk = Y{i};
            y_ijk(~obsIdx) = muSubj(~obsIdx) + (mosigmak * (oosigmak \ (yobs - muSubj(obsIdx))'))';

            % update mu_k. need to invert it back
            y_ijkr = g * y_ijk(:);
            Yhat(i, :, k) = y_ijkr';
        end
    end
    
    % mu update (exact same as in latentMissing)
    params.mu = zeros(K, dHigh);
    for k = 1:K
        gnk = wg.expect.gammank(:, k);
        params.mu(k, :) = sum(bsxfun(@times, gnk, Yhat(:, :, k))) ./ sum(gnk);
    end
    
    % update sigma using new mu
    % prepare sigma as K cells. This is faster than a 3-D matrix [nData-by-nData-by-K] 
    % since indexing is much faster this way.
    sigma = arrayfunc(@(s) zeros(dHigh, dHigh), 1:K);
    sigmaPrevCell = dimsplit(3, wg.params.sigma);
    for i = 1:N
        obsIdx = ydsmasks{i};
        r = R.data(R.idx{i}, :);
        g = G.data(:, G.idx{i});

        for k = 1:K
            gnk = wg.expect.gammank(i, k);
            % gamma rotation for missing values
            % gMissing rows won't add up to 1, but we only multiply it with a sparse S which would
            % have 0s in the gObserved entries.
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

            % empirical covariance for this point.
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
