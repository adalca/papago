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
    assert(dLow == wg.opts.model.dopca, ...
        'W low dim %d does not match explicit dopca parameter %d', dLow, wg.opts.model.dopca);

    % initialize new stats
    mu = zeros(K, dHigh);
    newW = zeros(dHigh, dLow, K);
    newsigmasq = zeros(1, K);
    
    doOld = false;
    if doOld
        muOld = zeros(K, dHigh);
        newWOld = zeros(dHigh, dLow, K);
    end
    
    
    % do everything *within* each cluster.
    for k = 1:K
        
        % TODO. testing various hacks.
        if isfield(data, 'estepW') && numel(wg.stats) < 10
            warning('mstep hack'); 
            obsMask = data.estepW == 1; 
        end
        
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
        
        if isfield(data, 'hackx') && data.hackx && numel(wg.stats) <= 1
            warning('HACKING X WITH RANDOM INIT');
            Xhat = randn(size(Xhat));
            Shat = repmat(eye(dLow), [1, 1, nEff]);
        end
        
        
        % TODO. testing various hacks.
        if isfield(data, 'estepW') && numel(wg.stats) < 10
            warning('mstep hack (end)'); 
            obsMask = data.W == 1; 
        end
        
        % compute old mu update (just for comparison)
        if doOld
            warning('checking old mu...');
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
                muOld(k, obsIdx) = muOld(k, obsIdx) + gnk .* muterm';

                % denominator for mu_k
                muDenom(obsIdx) = muDenom(obsIdx) + gnk;
            end
            muOld(k, :) = muOld(k, :) ./ muDenom;
        end
        
        % update W
        %   first, pre-compute (Shat * gammank), much faster than re-computing that operation for
        %   overlapping subsets of Shat inside the loop!
        gShatFull = bsxfun(@times, Shat, permute(gEff, [2,3,1]));
        hacnr = 0;
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
            wDenom = sgXhat' * sgXhat + sum(gShat, 3); % low-dim variance
            
            % hack add light random offsets to mu and W.
            if rcond(wDenom) < 1e-12
                hacnr = hacnr + 1;
                mu(k, j) = randn*0.001;
                newW(j, :, k) = randn([1, size(newW, 2)]) * 0.001;
                continue;
            end
            
            % sum over xbar
            xbar = xhat' * gnk;
            
            % update (closed-form) mu
            yBarJ = gnk' * yobs;
            yJsel = (gnk .* yobs)' * xhat * (wDenom \ xbar);
            bottom = sum(gnk) - gnk' * xhat * (wDenom \ xbar);
            top = (yBarJ - yJsel);
            mu(k, j) = top/bottom;
            
            % update W
            wNum = xhat' * (gnk .* (yobs - mu(k, j))); % low-dim centered data
            newW(j, :, k) = wNum' / wDenom;
            
            
            % old W
            if doOld
                wNum = xhat' * (gnk .* (yobs - muOld(k, j))); % low-dim centered data
                newWOld(j, :, k) = wNum' / wDenom;
            end
        end
        if hacnr > 0, fprintf('warning: hack rand added %d times due to bad rcond\n', hacnr); end
        assert(isclean(mu))
        if ~(isclean(newW))
            warning('mstepLS: found %d NANs in W_%d! Filling them in randomly', sum(sum(isnan(newW(:, :, k)))), k);
            newW(isnan(newW)) = randn(sum(isnan(newW(:))), 1);
        end
        
        % update residual variance
        v = 0;
        muDenom = zeros(1, dHigh);
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
            
            % updates
            v = v + resTerm * resTerm' + varTerm;
            muDenom(obsIdx) = muDenom(obsIdx) + gnk;
        end
        newsigmasq(k) = v ./ sum(muDenom);
        assert(isclean(newsigmasq), 'sigmasq is not clean');
        
    end
        
    % update parameters
    params = struct();
    params.mu = mu;
    params.W = newW;
    params.sigmasq = newsigmasq;

    
    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi), 'pi is not clean');
    