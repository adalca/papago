function params = mstepLatentSubspace(wg, data)

    % data
    Y = data.Y;
    if isfield(data, 'W')
        obsMask = data.W;
        assert(islogical(obsMask) | all(obsMask(:) == 0 | obsMask(:) == 1));
        obsMask = obsMask == 1; % in case it's not logical
    else
        obsMask = ~isnan(Y);
    end
    
    % dimensions and such
    N = size(data.Y, 1);
    K = data.K;
    Nk = sum(wg.expect.gammank, 1);
    [dHigh, dLow, ~] = size(wg.params.W);

    % initialize new stats
    mu = zeros(K, dHigh);
    newW = zeros(dHigh, dLow, K);
    newsigmasq = zeros(1, K);
    
    % do everything *within* each cluster.
    for k = 1:K

        % update Xhat and Shat in atlas space. 
        %   This should really be part of the E-step, but we'll do it here to
        %   not pass large matrices around. Specifically, here we only need to compute Shat and Xhat
        %   for one cluster at a time.
        Shat = zeros(dLow, dLow, N);
        Xhat = zeros(N, dLow);
        
        % go through each point
        for i = 1:N
            % don't need to compute Xhat_ki, Shat_ki when gnk is ~0
            % since they're only used when multiplied by gnk
            if wg.expect.gammank(i, k) < 1e-4, continue; end 
            
            % extract observed entries
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % prepare weights
            vk = wg.params.sigmasq(k);
            w = wg.params.W(obsIdx, :, k);
            L = w ./ sqrt(vk);

            % expectation updates
            %   ppca() uses Sherman-Morrison, probably because it's safer?
            %   S_ki = inv(L * L' + eye(dLow));
            ltl = L'*L;
            Ski = eye(dLow) - ltl / (eye(dLow) + ltl);
            Shat(:, :, i) = Ski;
            Xhat(i, :) = Ski/vk * (w' * (yobs - wg.params.mu(k, obsIdx))');
        end

        % update mu
        muDenom = zeros(1, dHigh);
        for i = 1:N
            % don't need to compute Xhat_ki, Shat_ki, since they're only used when gamma_ki is > 0
            gnk = wg.expect.gammank(i, k);
            if gnk < 1e-4, continue; end 
            
            % extract available voxels
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % extrace necessary data
            w = wg.params.W(obsIdx, :, k);
            X_ki = Xhat(i, :)';

            % update mu
            muterm = yobs' - w * X_ki;
            mu(k, obsIdx) = mu(k, obsIdx) + gnk .* muterm';

            % denominator for mu_k
            muDenom(obsIdx) = muDenom(obsIdx) + gnk;
        end
        % normalize mu
        % mutest = nanmean(Y' - wg.params.W * Xhat', 2);
        mu(k, :) = mu(k, :) ./ muDenom;
        if ~(isclean(mu))
            warning('mstepLS: found %d NANs in mu_%d! Filling them in', sum(isnan(mu(k, :))), k);
            allmu = mean(Y, 1);
            mu(k, isnan(mu(k, :))) = allmu(isnan(mu(k, :)));
        end

        % update W
        %   first, pre-compute (Shat * gammank), much faster than re-computing that operation for
        %   overlapping subsets of Shat inside the loop!
        gammank = permute(wg.expect.gammank(:, k), [2, 3, 1]);
        gShatFull = bsxfun(@times, Shat, gammank);
        for j = 1:dHigh
            % extract observed data indices
            % note: using gnk >= 1e-4
            lobsIdx = obsMask(:, j) & (wg.expect.gammank(:, k) >= 1e-4);

            % extract wanted data
            gnk = wg.expect.gammank(lobsIdx, k);
            xhat = Xhat(lobsIdx, :);
            gShat = gShatFull(:, :, lobsIdx); % Much faster this way: use pre-computed Shat*g
            
            % compute numerator and denominator
            sgXhat = bsxfun(@times, sqrt(gnk), xhat);
            wNum = xhat' * (gnk .* (Y(lobsIdx, j) - mu(k, j))); % low-dim centered data
            wDenom = sgXhat' * sgXhat + sum(gShat, 3); % low-dim variance
            
            % update W
            newW(j, :, k) = wDenom \ wNum;
        end
        if ~(isclean(newW))
            warning('mstepLS: found %d NANs in W_%d! Filling them in randomly', sum(sum(isnan(newW(:, :, k)))), k);
            newW(isnan(newW)) = randn(sum(isnan(newW(:))), 1);
        end
            
        % update residual variance
        v = zeros(dHigh, 1);
        for i = 1:N
            gnk = wg.expect.gammank(i, k);
            if gnk < 1e-4, continue; end
            
            % extract observed
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % new w
            w = newW(obsIdx, :, k);
            xhat = Xhat(i, :);
            
            % update terms
            resTerm = gnk .* (yobs - xhat * w' - mu(k, obsIdx)) .^ 2;
            varTerm = gnk .* diag(w * Shat(:, :, i) * w');

            v(obsIdx) = v(obsIdx) + resTerm(:) + varTerm(:);
        end
        newsigmasq(k) = sum(v) ./ sum(muDenom);
        assert(isclean(newsigmasq));
    end
        
    % update parameters
    params.mu = mu;
    params.W = newW;
    params.sigmasq = newsigmasq;
    sigma = zeros(dHigh, dHigh, K);
    for k = 1:K
        sigma(:,:,k) = newW(:,:,k) * newW(:,:,k)' + newsigmasq(k) .* eye(dHigh);
    end
    [params.sigma, params.sigmainv] = wgmm.sigmafull(wg, sigma);
    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));
    