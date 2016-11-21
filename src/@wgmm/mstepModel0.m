function params = mstepModel0(wg, data)

    % data
    wts = data.W;
    Y = data.Y;
    K = size(wg.expect.gammank, 2);
    N = size(Y, 1);
    Nk = sum(wg.expect.gammank, 1);
    assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
    
    mu = zeros(size(wgmm.params.mu));   
    % compute mu without weights
    for k = 1:K
        mu(k, :) = wg.expect.gammank(:, k)' * Y ./ sum(wg.expect.gammank(:, k));
    end
    
    % compute sigma without weights
    sigma = zeros(nDims, nDims, K);
    for k = 1:K
        xc = bsxfun(@minus, Y, mu(k, :));

        % Variant 1 -- this is slower
        % q = 0;
        % for i = 1:size(gammank, 1)
        %     q = q + gammank(i, k)' * (xc(:, i)' * xc(:, i));
        % end
        % sigma(:,:,k) = 1 ./ sum(gammank(:, k)) * q;

        % Variant 2 -- this is significantly faster!
        gw = sqrt(wg.expect.gammank(:, k));
        wx = bsxfun(@times, gw, xc);
        aa = wx' * wx;
        sigma(:, :, k) = aa ./ sum(wg.expect.gammank(:, k));
    end
    
    [params.sigma, params.sigmainv] = sigmafull(wg, sigma);
    
    params.mu = mu;

    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));