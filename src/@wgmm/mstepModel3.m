function params = mstepModel3(wg, data)
% latent missing data assumes just Y

    % data
    wts = data.W;
    Y = data.Y;
    K = size(wg.expect.gammank, 2);
    [N, dHigh] = size(Y);
    assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));

    mu = zeros(size(wgmm.params.mu));   
    for k = 1:K
        zz = bsxfun(@times, wg.expect.gammank(:, k), wts);
        mu(k, :) = sum(zz .* Y, 1) ./ sum(zz, 1);
    end

    % model 3, as computed in wgmm timeline. using matrix algebra, which seems significant
    % faster than the memsage method. It doesn't use too much memory, since we don't have to
    % compute as large matrices.
    wtw = W' * W;
    for k = 1:K
        xc = bsxfun(@minus, Y, mu(k, :));

        % Variant 1 -- this is slower
        % wx = W .* xc;
        % aa = bsxfun(@times, gammank(:, k), wx)' * wx;
        % wtsum = bsxfun(@times, gammank(:, k), W)' * W;
        % sigma(:, :, k) = aa ./ wtsum;

        % Variant 2 -- this is significantly faster!
        gw = bsxfun(@times, sqrt(wg.expect.gammank(:, k)), wts);
        wx = gw .* xc;
        aa = wx' * wx;
        wtsum = gw' * gw;
        sigma(:, :, k) = aa ./ wtsum;
        
        
        % reconstruction
        sigmarecon(:,:,k) = wgmm.sigmarecon(sigma(:,:,k), wtw, wg.opts.model.recon);

        % merge core with reconstruction
        switch wg.opts.model.merge
            case 'wfact'
                margs = {wg.opts.mergeargs};
            case 'wfact-mult'
                margs = {wg.opts.mergeargs};
            case 'wfact-mult-adapt'
                margs = {wg.opts.mergeargs{:}, size(X, 1), sum(gammank(:, k))};
            case 'freq-prior'                    
                margs = {wg.opts.mergeargs{:}};
            case 'none'
                margs = {};
            otherwise
                error('wgmm.sigmafull: Unknown combo method');
        end

        % merge sigmas
        sigmamerge(:,:,k) = wgmm.sigmamerge(sigma(:,:,k), sigmarecon(:,:,k), ...
            wtw, wg.opts.model.merge, margs{:});
        sigma(:,:,k) = sigmamerge(:,:,k);
    end          
    
    

    
    
    params.mu = mu;
    [params.sigma, params.sigmainv] = sigmafull(wg, sigma);
    
    
    % 'memsafe-model3'
    % model3, with internal loop (memory safe). 
    % Can also be used to make sure model1 implementation above yields similar results.
%     for k = 1:K
%         numer = 0;
%         denom = 0;
%         for i = 1:size(X, 1)
%             w = W(i, :)';
%             gsm = gammank(i, k) * w;
%             numer = numer + gsm .* X(i, :)';
%             denom = denom + gsm;
%         end
%         mu(k, :) = numer ./ denom;
%     end
%
%     % memory safe (looping) model 3
%     for k = 1:K
%         xc = bsxfun(@minus, X, mu(k, :))';
%         numer = 0;
%         denom = 0;
%         for i = 1:size(X, 1)
%             w = sqrt(W(i, :)');
%             z = w .* xc(:, i);
%             q = z * z';
%             numer = numer + gammank(i, k) * q;
% 
%             denom = denom + gammank(i, k) .* (w * w');
%         end
%         sigma(:, :, k) = numer ./ denom;
%     end
