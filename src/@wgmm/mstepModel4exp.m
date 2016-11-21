function params = mstepModel4exp(wg, data)
% latent missing data assumes just Y

    % data
    wts = data.W;
    Y = data.Y;
    K = size(wg.expect.gammank, 2);
    [N, dHigh] = size(Y);
    assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
    
    % mu update
    mu = zeros(size(wgmm.params.mu));   
    for k = 1:K
        numer = 0;
        denom = 0;
        sigmak = wgmm.params.sigma(:,:,k);
        for i = 1:size(Y, 1)
            w = W(i, :);
            yobs = Y(i, :);

            % sigmas
            Di = wgmm.model4fn(w);
            sigma = sigmak + Di;

            numer = numer + gammank(i, k) * ((sigma + Di) \ yobs(:));
            denom = denom + gammank(i, k) * inv(sigma + Di);
        end
        mu(k, :) = denom \ numer; 
    end

    % sigma update
    for k = 1:K
        xc = bsxfun(@minus, Y, mu(k, :))';

        numer = 0;
        denom = 0;
        for i = 1:N
            w = sqrt(W(i, :)');
            Di = modelopts.model4fn(w);
            s = wg.sigma(:,:,k);
            sg = s + Di;
            xtx = xc(:, i) * xc(:, i)';

            numer = numer + gammank(i, k) * (xtx + (s / sg) * Di);

            denom = denom + gammank(i, k);
        end
        sigma(:, :, k) = numer ./ denom;
        imagesc(sigma(:, :, k)); 
        drawnow;
        pause(0.01);
    end