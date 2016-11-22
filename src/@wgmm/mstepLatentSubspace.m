function params = mstepLatentSubspace(wg, data)

    % data
    wts = data.W;
    Y = data.Y;
    
    % dimensions and such
    K = data.K;
    Nk = sum(wg.expect.gammank, 1);
    [dHigh, dLow, ~] = size(wg.params.W);

    % update Xhat in atlas space. 
    %   This should really be part of the E-step, but we'll do it here to
    %   not pass large matrices around. Maybe we should fix this for clenliness!
    muk = zeros(K, size(Y,2));
    denom = zeros(K, size(Y,2));
    Shat = zeros(dLow, dLow, size(Y,1), K);
    Xhat = zeros(size(Y,1), dLow, K);
    for i = 1:size(Y, 1)
        obsIdx = wts(i, :) == 1;
        yobs = Y(i, obsIdx);

        for k = 1:K
            sigmasq_t = wg.params.sigmasq(k);
            w = wg.params.W(obsIdx, :, k);
            L = w ./ sqrt(sigmasq_t);

            % expectation updates. We cannot store these reliable
            % since they would be very large.
            % ppca() uses Sherman-Morrison, probably because it's safer?
            % S_ki = inv(L * L' + eye(dLow));
            ltl = L'*L;
            S_ki = eye(dLow) - ltl / (eye(dLow) + ltl);
            X_ki = S_ki/sigmasq_t * (w' * (yobs - wg.params.mu(k, obsIdx))');

            % update mu_k
            muterm = yobs' - w * X_ki;
            muk(k, obsIdx) = muk(k, obsIdx) + wg.expect.gammank(i, k) .* muterm';

            % denominator for mu_k
            denom(k, obsIdx) = denom(k, obsIdx) + wg.expect.gammank(i, k);

            % update large matrices
            Shat(:,:,i,k) = S_ki;
            Xhat(i,:,k) = X_ki';
        end
    end

    % normalize mu
    % mutest = nanmean(Y' - wg.params.W * Xhat', 2);
    mu = muk./denom;

    % update W
    % have to recompute S_ki and X_ki again. This is slow. Maybe
    % try to store them in memory :(
    newW = zeros(dHigh, dLow, K);
    for j = 1:dHigh
        obsIdx = wts(:, j) == 1;
        yobs = Y(obsIdx, j);

        for k = 1:K
            x = Xhat(obsIdx, :, k);
            g = wg.expect.gammank(obsIdx, k);

            gx = bsxfun(@times, g, x);
            gS = bsxfun(@times, Shat(:, :, obsIdx, k), permute(g, [2, 3, 1]));
            t1 = gx' * gx + sum(gS, 3);

            t2 = x' * (g .* (yobs - mu(k, j)));

            newW(j, :, k) = t1 \ t2;
        end
    end

    % update residual variance
    t1a = zeros(dHigh, K);
    t2a = zeros(dHigh, K);
    for i = 1:size(Y, 1)
        obsIdx = wts(i, :) == 1;
        yobs = Y(i, obsIdx);

        for k = 1:K
            w = newW(obsIdx, :, k);
            x = Xhat(i, :, k);
            g = wg.expect.gammank(i, k);

            t1 = g .* (yobs - x * w' - mu(k, obsIdx)) .^ 2;
            t2 = g .* diag(w * Shat(:, :, i, k) * w');

            t1a(obsIdx, k) = t1a(obsIdx, k) + t1(:);
            t2a(obsIdx, k) = t2a(obsIdx, k) + t2(:);
            % v(obsIdx, k) = v(obsIdx, k) + t1(:) + t2(:);
        end
    end
    v = t1a + t2a;
    newsigmasq = sum(v, 1) ./ sum(denom, 2)';

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
    