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
    
    % prepare new variables
    newMu = zeros(K, dHigh);
    newSigma = zeros(dHigh, dHigh, K);
    
    % Operate per cluster since there is no interaction among clusters, and this way we can keep
    % much smaller variables around.
    for k = 1:K
        sigmaPrev = wg.params.sigma(:, :, k);
        muPrev = wg.params.mu(k, :);
        expectk = wg.expect.gammank(:, k);

        % update Yhat 
        %   This should really be part of the E-step, but we'll do it here to
        %   not pass large matrices around. 
        %   y_ik^Oi (observed) = y_i^Oi data
        %   y_ik^Mi (missing) = mu_k^Mi + C_k^MiOi * inv(C_k^OiOi) * (y_i^Oi - mu_k^Oi)
        Yhat = Y; % N-by-dHigh;
        for i = 1:N
            % don't need to compute Yhat_ki, since it's only used when gamma_ki is > 0
            if expectk(i) < 1e-4, continue; end 
            
            % extract the observed y values
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % extract sigma
            oosigmak = sigmaPrev(obsIdx, obsIdx);
            mosigmak = sigmaPrev(~obsIdx, obsIdx);

            % update Yhat
            Yhat(i, ~obsIdx) = muPrev(~obsIdx) + (mosigmak * (oosigmak \ (yobs - muPrev(obsIdx))'))';
        end

        % update mu_k
        newMuk = sum(bsxfun(@times, expectk, Yhat)) ./ sum(expectk);

        % update sigma. Start with the empirical covariance. Note:
        %   - ECMNMLE use the most recent mu to *recompute* y_ik. This does not make sense to us
        %     y_ik = Y(i, :);
        %     y_ik(~obsIdx) = newMuk(~obsIdx) + (mosigmak * (oosigmak \ (Y(i, obsIdx) - newMuk(obsIdx))'))';
        ydiff = bsxfun(@times, sqrt(wg.expect.gammank(:, k)), bsxfun(@minus, Yhat, newMuk));
        newSigmak = ydiff' * ydiff;
        
        % add the correction term
        for i = 1:N
            if expectk(i) < 1e-4, continue; end
            
            % extract the observed mask
            obsIdx = obsMask(i, :);

            % extract appropriate sigmas
            oosigmak = sigmaPrev(obsIdx, obsIdx);
            mosigmak = sigmaPrev(~obsIdx, obsIdx);
            
            % compute correction term
            sCorrTerm1 = sigmaPrev(~obsIdx, ~obsIdx);
            sCorrTerm2 = mosigmak * (oosigmak \ mosigmak');
            sCorr = sCorrTerm1 - sCorrTerm2;

            % update sigma
            newSigmak(~obsIdx, ~obsIdx) = newSigmak(~obsIdx, ~obsIdx) + expectk(i) * sCorr;
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
