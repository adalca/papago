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
    assert(dLow == wg.opts.model.dopca);

    % initialize new stats
    mu = zeros(K, dHigh);
    newW = zeros(dHigh, dLow, K);
    newsigmasq = zeros(1, K);
    
    % do everything *within* each cluster.
    for k = 1:K
        
        % compute which gammank are non-trivial for this cluster. For many operations, we only need
        % to work for those elements where gamma_nk > 0, saving a good amount of computation and,
        % perhaps more importantly, memory.
        gEffLog = wg.expect.gammank(:, k) >= 1e-4;
        gEffIdx = find(gEffLog);
        gEff = wg.expect.gammank(gEffIdx, k);
        nEff = numel(gEffIdx);

        % update Xhat and Shat in atlas space. 
        %   This should really be part of the E-step, but we'll do it here to
        %   not pass large matrices around. Specifically, here we only need to compute Shat and Xhat
        %   for one cluster at a time.
        Shat = zeros(dLow, dLow, nEff);
        Xhat = zeros(N, dLow);
               
        for ei = 1:nEff
            i = gEffIdx(ei);
            
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
            Shat(:, :, ei) = Ski;
            Xhat(ei, :) = Ski/vk * (w' * (yobs - wg.params.mu(k, obsIdx))');
        end
        assert(isclean(Xhat), 'Xhat is not clean :(');
        assert(isclean(Shat), 'Shat is not clean :(');

        % update mu
        muDenom = zeros(1, dHigh);
        for ei = 1:nEff
            i = gEffIdx(ei);
            gnk = gEff(ei);
            
            % extract available voxels
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % extrace necessary data
            w = wg.params.W(obsIdx, :, k);
            X_ki = Xhat(ei, :)';

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
        gShatFull = bsxfun(@times, Shat, permute(gEff, [2,3,1]));
        for j = 1:dHigh
            % extract observed data indices
            lobsIdx = obsMask(gEffIdx, j); 
            yobs = Y(gEffIdx(lobsIdx), j);

            % extract wanted data
            gnk = gEff(lobsIdx);
            xhat = Xhat(lobsIdx, :);
            gShat = gShatFull(:, :, lobsIdx); % Much faster this way: use pre-computed Shat*g
            
            % compute numerator and denominator
            sgXhat = bsxfun(@times, sqrt(gnk), xhat);
            wNum = xhat' * (gnk .* (yobs - mu(k, j))); % low-dim centered data
            wDenom = sgXhat' * sgXhat + sum(gShat, 3); % low-dim variance
            
            % update W
            newW(j, :, k) = wDenom \ wNum;
        end
        if ~(isclean(newW))
            warning('mstepLS: found %d NANs in W_%d! Filling them in randomly', sum(sum(isnan(newW(:, :, k)))), k);
            newW(isnan(newW)) = randn(sum(isnan(newW(:))), 1);
        end
        
        % update residual variance
        v = 0;
        for ei = 1:nEff
            i = gEffIdx(ei);
            gnk = gEff(ei);
            
            % extract observed
            obsIdx = obsMask(i, :);
            yobs = Y(i, obsIdx);

            % new w
            w = newW(obsIdx, :, k);
            xhat = Xhat(ei, :);
            
            % update terms
            resTerm = sqrt(gnk) .* (yobs - xhat * w' - mu(k, obsIdx)); % will square below
            % want to compute the term sum(gnk .* diag(w * Shat(:, :, ei) * w');)
            % which is really the trace, so we can change out multiplications around to be faster.
            m1 = w' * w * Shat(:, :, ei); % fastest version for trace.
            varTerm = gnk .* trace(m1);
            
            v = v + resTerm * resTerm' + varTerm;
        end
        newsigmasq(k) = v ./ sum(muDenom);
        assert(isclean(newsigmasq));
        
    end
        
    % update parameters
    params = struct();
    params.mu = mu;
    params.W = newW;
    params.sigmasq = newsigmasq;

    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));
    