function wg = init(wg, X, W, K, varargin)
% initializing a wgmm fit optimization. 

    assert(size(X, 1) >= K);
    assert(all(size(X) == size(W)));
    
    % prepare convenient variables
    [N, D] = size(X);   

    switch wg.initmethod
        
        % initialize from another wgmm
        case 'wgmm' 
            wgm = varargin{1};
            wg.mu = wgm.mu;
            wg.sigma = wgm.sigma;
            wg.pi = wgm.pi;
            
        % initialize from a matlab gmm
        case 'gmm'
            gm = varargin{1};
            wg.mu = gm.mu;
            wg.sigma = gm.Sigma;
            wg.pi = gm.ComponentProportion;
        
        % assign a cluster number to each point randomly
        case 'randset' 
            
            initv = zeros(1, N);
            
            % make sure there is at least one element in each cluster
            initv(1:K) = 1:K; 
            initv = initv(randperm(N));
            
            % assign clusters
            initv(initv == 0) = randi([1, K], 1, N - K);
            
            % compute initial parameters
            for k = 1:K
                x = X(initv == k, :);
                wg.mu(k, :) = mean(x, 1);
                wg.sigma(:,:,k) = cov(x);
                wg.sigmainv(:,:,k) = inv(wg.sigma(:,:,k));
                wg.pi(k) = numel(x) ./ N;
            end
            
        % initialize exemplars as means. Sigma is a diagnoal of variances, and pi is 1/K
        case 'exemplar' 
            p = randsample(N, K);
            
            % compute initial parameters
            for k = 1:K
                wg.mu(k, :) = X(p(k), :);
                
                wg.sigma(:, :, k) = diag(var(X) + wg.sigmareg);
                wg.sigmainv(:, :, k) = inv(wg.sigma(:,:,k));
                wg.pi(k) = 1 ./ K;
            end
            
        % initialize exemplars as means. 
        % sigma and pi are learned via cluster assignments based on the 'exemplar' initialization
        case 'exemplar-spec' 
            p = randperm(N);
            
            
            % compute initial parameters
            for k = 1:K
                wg.mu(k, :) = X(p(k), :);
                
                wg.sigma(:, :, k) = diag(var(X));
                wg.sigmainv(:, :, k) = inv(wg.sigma(:,:,k));
                wg.pi(k) = K ./ N;
            end
            
            n = wg.logpost();
            [~, mi] = min(n, [], 2);
            for k = 1:K
                x = X(mi == k, :);
                
                params.sigma(:, :, k) = cov(x);
                params.sigmainv(:, :, k) = inv(params.sigma(:,:,k));
                params.pi(k) = sum(mi == k) ./ N;
            end
            
        case 'git-gmm'

            wgm = fitgmdist(X, K, varargin{:});
            
            wg.mu = wgm.mu;
            wg.pi = wgm.ComponentProportion;
            wg.sigma = wgm.Sigma;
            for k = 1:K
                wg.sigmainv(:,:,k) = inv(wg.sigma(:,:,k));
            end
            
        case 'model4exp05'
            assert(K == 1)
            wg.mu = nanmean(X);
            for i = 1:size(X, 1)
                X(i, W(i, :) < 0.5) = wg.mu(W(i,:) < 0.5);
            end
            wg.sigma = cov(X);
            wg.pi = 1;
            
            
        otherwise
            error('wgmm.init: unknown init method');
    end
    
    % check cleanliness and sizes.
    assert(isclean(wg.mu));
    assert(isclean(wg.sigma));
    assert(isclean(wg.pi));
    assert(all(size(wg.mu) == [K, D]), 'The size of init mu is incorrect');
    assert(all(size(permute(wg.sigma, [3 1 2])) == [K, D, D]), 'The size of init sigma is incorrect');
    assert(numel(wg.pi) == K, 'The size of init pi is incorrect');
end
