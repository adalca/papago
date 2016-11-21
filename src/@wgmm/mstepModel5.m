function params = mstepModel5(wg, data)
% latent missing data assumes just Y

    % data
    wts = data.W;
    Y = data.Y;
    K = size(wg.expect.gammank, 2);
    [N, dHigh] = size(Y);
    assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
    
    % E-step. Should compute this in Estep.
    murnk = zeros(N, K, D); % will get permuted to N-D-K
    for k = 1:K
        muk = wgmm.params.mu(k, :);
        sigmak = wgmm.params.sigma(:,:,k);

        for i = 1:N
            w = W(i, :);
            x = X(i, :);

            % sigmas
            Di = wgmm.opts.model.model4fn(w);
            sigma = sigmak + Di;

            % z = (sigmak / sigma) * (x - muk)' + muk';
            z = sigmak * (sigma \ (x - muk)') + muk'; % faster
            murnk(i, k, :) = reshape(z, [1, 1, D]);
        end 
    end
    assert(isclean(murnk), 'murnk is unclean');
    
    % compute mu without weights
    Xr = permute(murnk, [1, 3, 2]); % is N-by-D-by-K. Added this much later. unsure if correct. It's a mess, basically.
    mu = zeros(size(wgmm.params.mu));   
    for k = 1:K
        mu(k, :) = gammank(:, k)' * Xr(:,:,k) ./ sum(gammank(:, k));
    end
    

    %         for k = 1:K
    %             muk = wg.mu(k, :);
    %             sigma(:,:,k) = fminsearch(@(s) wg.model5exp(s, X, W, muk, wg.model4fn), wg.sigma(:,:,k));
    %         end
    % 
    %         for k = 1:K
    %             numer = 0;
    %             denom = 0;
    %             sigmak = wg.sigma(:,:,k);
    %             for i = 1:size(X, 1)
    %                 w = W(i, :);
    %                 x = X(i, :);
    %                 df = (x + wg.mu(k,:));
    %                 dfd = df(:) * (df(:))';
    % 
    %                 % sigmas
    %                 Di = wg.model4fn(w);
    %                 sg = sigmak + Di;
    % 
    %                 numer = numer + gammank(i, k) * dfd / (sg + Di) / Di / (inv(sigmak) + inv(Di));
    %                 denom = denom + gammank(i, k) * inv(sg + Di);
    %             end
    %             sigma(:, :, k) = denom \ numer; 
    %         end    
    % 
    %         for k = 1:K
    %             numer = 0;
    %             denom = 0;
    %             sigmak = wg.sigma(:,:,k);
    %             for i = 1:size(X, 1)
    %                 w = W(i, :);
    %                 x = X(i, :);
    %                 df = (x - wg.mu(k,:));
    %                 dfd = df(:) * (df(:))';
    % 
    %                 % sigmas
    %                 Di = wg.model4fn(w);
    %                 sg = sigmak + Di;
    % 
    %                 numer = numer + gammank(i, k) * (- eye(numel(x)) + (dfd / (sg + Di))) / (sg + Di);
    %             end
    %             sigma(:, :, k) = sigmak - 0.00001 * numer ./ sum(gammank(:, k)); 
    %         end    
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % MICCAI2016 Push implementation.
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     fprintf(2, 'current using init with sigma of y (noisy data). Maybe use previous sigma?');
    % 
    %     for k = 1:K
    %     muk = wg.mu(k, :);
    %     xc = bsxfun(@minus, X, muk);
    % 
    %     % compute sigma for this cluster ffrom noisy data xc.
    %     gw = sqrt(gammank(:, k));
    %     wx = bsxfun(@times, gw, xc);
    %     aa = wx' * wx; % gamma * (x-u) (x-u)'
    %     sigmainitk = aa ./ sum(gammank(:, k));
    % 
    %     % hack 1
    %     numer = 0; denom = 0; 
    %     meanD = 0;
    %     for i = 1:size(X, 1)
    %         w = W(i, :); 
    %         x = X(i, :);
    %         df = (x - muk);
    %         dfd = df(:) * (df(:))';
    % 
    %         Di = wg.model4fn(w);
    %         sg = sigmainitk + Di;
    % 
    %         numer = numer + (sg \ dfd);
    %         %denom = denom + inv(sg);
    %         denom = denom + inv(sg) + inv(sg) * dfd * inv(sigmainitk) * Di * inv(sigmainitk);
    % 
    %         meanD = meanD + Di;
    %     end
    %     meanD = meanD ./ size(X, 1);
    %     sigma(:,:,k) = (denom \ numer) - meanD ./ 8; % empirical ./4
    %     % sigma(:,:,k) = (denom \ numer); 

    % new implementation
    for k = 1:K
        muk = wg.mu(k, :);

        % model 5
        df = bsxfun(@minus, wg.expect.Xk(:,:,k), muk);
        numer = 0; 
        for i = 1:size(X, 1)
            % extract important terms.
            w = W(i, :);
            Di = modelopts.model4fn(w);
            dfd = df(i, :)' * df(i, :);
            sk = wg.sigma(:,:,k);

            % previous sigma^{Rt}_{ik}. We can't precompute this due to the size.
            sigmaikr = Di * ( (Di + sk) \ sk);

            % update numerator and denominator.
            numer = numer + gammank(i, k) .* (sigmaikr + dfd);
        end
        sigma(:,:,k) = numer ./ sum(gammank(:, k));
    end
    
    params.mu = mu;
    params.sigma = sigma;
    
    