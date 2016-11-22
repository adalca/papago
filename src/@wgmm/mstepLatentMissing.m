function params = mstepLatentMissing(wg, data)
% latent missing data assumes just Y

    % data
    Y = data.Y;
    if isfield(data, 'W')
        obsMask = data.W;
        assert(islogical(obsMask) | all(obsMask(:) == 0 | obsMask(:) == 1));
        obsMask = obsMask == 1; % in case it's not logical
    else
        obsMask = ~isnan(Y);
    end
    
    % dimensions and numbers
    K = data.K;         
    [N, dHigh] = size(Y);
    Nk = sum(wg.expect.gammank, 1);
    
    % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than indexing
    % into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times per iteration, it
    % ends up mattering for our iterations.
    sigmaPrevCell = dimsplit(3, wg.params.sigma);
    
    % update Yhat 
    %   This should really be part of the E-step, but we'll do it here to
    %   not pass large matrices around. Maybe we should fix this for clenliness!
    %   y_ik^Oi (observed) = y_i^Oi data
    %   y_ik^Mi (missing) = mu_k^Mi + C_k^MiOi * inv(C_k^OiOi) * (y_i^Oi - mu_k^Oi)
    Yhat = zeros(N, dHigh, K);
    for i = 1:N
        % extract the observed y values
        obsIdx = obsMask(i, :);
        yobs = Y(i, obsIdx);
        
        % go through each cluster
        for k = 1:K
            % extract sigma and mu
            sigmaPrev = sigmaPrevCell{k};
            oosigmak = sigmaPrev(obsIdx, obsIdx);
            mosigmak = sigmaPrev(~obsIdx, obsIdx);
            mu = wg.params.mu(k, :);

            % update Y
            ymis = mu(~obsIdx) + (mosigmak * (oosigmak \ (yobs - mu(obsIdx))'))';
            Yhat(i, ~obsIdx, k) = ymis;
        end
    end

    % update mu
    params.mu = zeros(K, dHigh);
    for k = 1:K
        gnk = wg.expect.gammank(:, k);
        params.mu(k, :) = sum(bsxfun(@times, gnk, Yhat(:, :, k))) ./ Nk(k);
    end
    
    % update sigma using new mu
    sigma = arrayfunc(@(s) zeros(dHigh, dHigh), 1:K);
    for i = 1:N
        % extract the observed mask
        obsIdx = obsMask(i, :);
        
        % go through each cluster
        for k = 1:K
            gnk = wg.expect.gammank(i, k);
            
            % compute S correction term
            sigmaPrev = sigmaPrevCell{k};
            % extract appropriate sigmas
            oosigmak = sigmaPrev(obsIdx, obsIdx);
            mosigmak = sigmaPrev(~obsIdx, obsIdx);
            % compute correction term
            sCorrTerm1 = sigmaPrev(~obsIdx, ~obsIdx);
            sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
            sCorr = sCorrTerm1 - sCorrTerm2;

            % prep empirical covariance
            y_ik = Yhat(i, :, k);
            sEmpirical = (y_ik - params.mu(k, :))' * (y_ik - params.mu(k, :));

            % update sigma
            sigma{k}(~obsIdx, ~obsIdx) = sigma{k}(~obsIdx, ~obsIdx) + gnk * sCorr;
            sigma{k} = sigma{k} + gnk * sEmpirical;
        end
    end
    params.sigma = cat(3, sigma{:});
    assert(isclean(params.sigma), 'sigma is not clean');

    % normalize
    for k = 1:K
        params.sigma(:, :, k) = params.sigma(:, :, k) ./ Nk(k);
    end

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
            W = u(:, 1:dLow) * sqrt(s(1:dLow, 1:dLow) - sigmasq * eye(dLow));
            Covar = W*W' + sigmasq*eye(size(u,2));

            params.sigmasq(k) = sigmasq;
            params.W(:,:,k) = W;
            params.sigma(:,:,k) = Covar;
        end
    end
    
    [params.sigma, params.sigmainv] = wgmm.sigmafull(wg, params.sigma);
    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));

%     else % use low-dimensional model for y
%         % this assumes y = W'x + mu + eps, for some low-dimensional x (that is not
%         % explicitly modelled). This does not change the update, but C is W'W + sI
%         % for W beind some [small x large] matrix, meaning that we can save on
%         % computation in updating y_ik^Mi if we are careful.
%         % y_ik^Mi = mu_k^Mi + C_k^MiOi * inv(C_k^OiOi) * (y_i^Oi - mu_k^Oi)
%         % ...
%         %    = mu_k^Mi + [Lm'Lo - Lm'LoLo' * inv(I+LoLo') * Lo] * (y_i^Oi - mu_k^Oi)
%         % where Lm are the missing columns of s*W, and Lo are the observed colums 
%         %
%         % This, in part, is using the trick of computing inverses using generalized
%         % sherman-morrison (formula in matrix cookbook 166) 
%         % inv(I + L'L) = I + L' inv(I+LL') L
%         % This is very useful is p << d.
% 
%         error('The pca hack is not much faster, but currently more prone to error');
% 
%         % get some data ready
%         assert(wgmm.opts.usepca == size(wgmm.params.W, 1));
%         L = 1./sqrt(wgmm.params.sigmasq) * wgmm.params.W; 
%         pcak = wgmm.opts.usepca;
% 
%         % compute
%         Lm = L(:, ~obsIdx);
%         Lo = L(:, obsIdx);
%         LooT = Lo * Lo';
%         t1 = Lm' * Lo - Lm' * LooT * (eye(pcak) + LooT) \ Lo;
%     end
