function params = mstepModel1(wg, data)

    % data
    wts = data.W;
    Y = data.Y;
    K = size(wg.expect.gammank, 2);
    Nk = sum(wg.expect.gammank, 1);
    gammank = wg.expect.gammank;
    assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
    
    % mu
    mu = zeros(size(wgmm.params.mu));   
    for k = 1:K
        % method 1. using matrix algebra, but running into memory issues, and not that much faster.
        wwtsigmainv = wgmm.iwAiw(wts, wgmm.params.sigmainv(:,:,k)); % warning: this is still too big.
        s = wgmm.sx(Y, wwtsigmainv);
        bottom = sum(bsxfun(@times, gammank(:,k), permute(wwtsigmainv, [3, 1, 2])), 1);
        mu = sum(bsxfun(@times, gammank(:, k), s), 1) / squeeze(bottom);
    end
    
    % model 1 in terms of derivations during the wgmm project. This implementation is using
    % matrix algebra, but running into memory issues, and not that much faster. see
    % memsafe-model1 method. This version uses Paolo de Leva's multiprod() function.
    iwAiw = coreargs.iwAiw;
    for k = 1:K
        xc = bsxfun(@minus, Y, mu(k, :))';

         % method 1. 
        diff = permute(xc, [1, 3, 2]);
        xxt = multiprod(diff, permute(diff, [2, 1, 3])); % maybe just loop?
        ixxti = iwAiw(wts, xxt); 
        s1 = sum(bsxfun(@times, ixxti, permute(gammank(:, k), [2, 3, 1])), 3);
        sigma = 1 ./ Nk(k) .* s1;
    end
    
    params.mu = mu;
    [params.sigma, params.sigmainv] = sigmafull(wg, sigma);

    % pi update
    params.pi = Nk ./ N; 
    assert(isclean(params.pi));
    
%     % model1, with internal loop (memory safe). 
%     % Can also be used to make sure model1 implementation above yields similar results.
%     for k = 1:K
%         numer = 0; 
%         denom = 0;
%         for i = 1:size(X, 1)
%             localwt = W(i, :)';
%             wwt = localwt * localwt';
%             sm = wwt .* wgmm.params.sigmainv(:,:,k);
%             gsm = gammank(i, k) * sm;
%             numer = numer + gsm * X(i, :)';
%             denom = denom + gsm;
%         end
%         mu(k, :) = denom \ numer;
%     end
%
%     % result should be the same as 'model1'
%     for k = 1:K
%         xc = bsxfun(@minus, X, mu(k, :))';
% 
%         % method 2. loop.
%         numer = 0; 
%         for i = 1:size(X, 1)
%             w = W(i, :)';
%             wx = w .* xc(:, i); 
%             q = wx * wx'; 
%             numer = numer + gammank(i, k) * q;
%         end
%         sigma(:,:,k) = numer ./ Nk(k);
%     end