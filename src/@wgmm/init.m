function wg = init(wg, data, varargin)
% initializing a wgmm fit optimization. 
    Y = data.Y;
    K = data.K;

    assert(size(Y, 1) >= K);
    
    if isfield(data, 'W')
        W = data.W;
        assert(all(size(Y) == size(W)));
    end
    
    % prepare convenient variables
    [N, D] = size(Y);   

    switch wg.opts.init.method
        
        % initialize from another wgmm
        case {'wgmm', 'wgmmr', 'wgmmI'}
            wg.params = varargin{1}.wgmm.params;
            % wg.expect = varargin{1}.wgmm.expect;
            
            if strcmp(wg.opts.init.method, 'wgmmr')
                wg.params.W = wg.params.W + randn(size(wg.params.W)) * 0.05;
                wg.params.sigma = wg.wv2sigma();
            end
            
            if strcmp(wg.opts.init.method, 'wgmmI')
                e = repmat(eye(D), [1,1,K]);
                wg.params.W = wg.params.W + e(:, 1:size(wg.params.W, 2), :);
                wg.params.sigma = wg.wv2sigma();
            end            
            
            if strcmp('latentSubspace', wg.opts.model.name) && ...
                (~isfield(wg.params, 'W') || size(wg.params.W, 2) ~= wg.opts.model.dopca)
                warning('wgmm init: W not the right size. Re-estimating from Sigma');

                W = zeros(size(Y, 2), wg.opts.model.dopca, K);
                v = zeros(1, K);
                for k = 1:K
                    [W(:,:,k), v(k)] = sigma2modelParams(wg.params.sigma(:,:,k), wg.opts.model.dopca);
                    wg.params.sigma(:,:,k) = W(:,:,k) * W(:,:,k)' + v(k) * eye(size(Y, 2));
                end
                wg.params.W = W;
                wg.params.sigmasq = v;
            end
            
            
            
        % initialize from a matlab gmm
        case 'gmm'
            gm = varargin{1};
            wg.params.mu = gm.mu;
            wg.params.sigma = gm.Sigma;
            wg.params.pi = gm.ComponentProportion;
        
        % assign a random cluster idx to each data point 
        case 'randset' 
            
            initv = zeros(1, N);
            
            % make sure there is at least one element in each cluster
            initv(1:K) = 1:K; 
            initv = initv(randperm(N));
            
            % assign clusters
            initv(initv == 0) = randi([1, K], 1, N - K);
            
            % compute initial parameters
            for k = 1:K
                x = Y(initv == k, :);
                w = W(initv == k, :);
                wg.params.mu(k, :) = wmean(x, w, 1);
                wg.params.sigma(:,:,k) = nancov(x);
                wg.params.sigmainv(:,:,k) = inv(wg.params.sigma(:,:,k));
                wg.params.pi(k) = numel(x) ./ N;
            end
            
        % initialize exemplars as means. Sigma is a diagnoal of variances, and pi is 1/K
        case 'exemplar' 
            p = randsample(N, K);
            
            % compute initial parameters
            for k = 1:K
                wg.params.mu(k, :) = Y(p(k), :);
                
                wg.params.sigma(:, :, k) = diag(var(Y) + wg.opts.regularizationValue);
                wg.params.sigmainv(:, :, k) = inv(wg.params.sigma(:,:,k));
                wg.params.pi(k) = 1 ./ K;
            end
            
        % initialize exemplars as means. 
        % sigma and pi are learned via cluster assignments based on the 'exemplar' initialization
        case 'exemplar-spec' 
            p = randperm(N);
            
            
            % compute initial parameters
            for k = 1:K
                wg.params.mu(k, :) = Y(p(k), :);
                
                wg.params.sigma(:, :, k) = diag(var(Y));
                wg.params.sigmainv(:, :, k) = inv(wg.params.sigma(:,:,k));
                wg.params.pi(k) = K ./ N;
            end
            
            n = wg.logpost();
            [~, mi] = min(n, [], 2);
            for k = 1:K
                x = Y(mi == k, :);
                
                params.sigma(:, :, k) = cov(x);
                params.sigmainv(:, :, k) = inv(params.sigma(:,:,k));
                params.pi(k) = sum(mi == k) ./ N;
            end
            
        case 'git-gmm'

            wgm = fitgmdist(Y, K, varargin{:});
            
            wg.params.mu = wgm.mu;
            wg.params.pi = wgm.ComponentProportion;
            wg.params.sigma = wgm.Sigma;
            for k = 1:K
                wg.params.sigmainv(:,:,k) = inv(wg.params.sigma(:,:,k));
            end
            
        case 'model4exp05'
            assert(K == 1)
            wg.params.mu = nanmean(Y);
            for i = 1:size(Y, 1)
                Y(i, W(i, :) < 0.5) = wg.params.mu(W(i,:) < 0.5);
            end
            wg.params.sigma = cov(Y);
            wg.params.pi = 1;
            
        case 'latentMissing-randIdx-twostage'
            % randomly choose cluster assignments and initate clusters 
            % in 'twostage' - Estimate mean, fill NaNs with mean, then estimate covar.
            
            % prepare useful variables
            obsIdx = W == 1;
            yfill = Y;
            yfill(~obsIdx) = nan;
            
            % choose random cluster assignments
            ridx = randi([1, K], [N, 1]);
            
            % go through each cluster and initialize in twostage
            sigmas = zeros(D, D, K);
            mu = zeros(K, D);
            for k = 1:K
                % find mean and fill in means
                mu(k,:) = nanmean(Y(ridx == k, :));
                for i = find(ridx == k)'
                    yfill(i, ~obsIdx(i,:)) = mu(k, ~obsIdx(i,:));
                    
                    % if 
                end
                
                % compute covariance
                sigmas(:,:,k) = cov(yfill(ridx == k, :));
            end
            
            % prepare init structure
            wg.params.mu = mu;
            wg.params.sigma = sigmas;
            wg.params.pi = hist(ridx, 1:K) ./ N;
            

            
        case 'latentMissing-clusterIdx-twostage'
            % given cluster assignments and initate clusters 
            % in 'twostage' - Estimate mean, fill NaNs with mean, then estimate covar.
            
            % prepare useful variables
            obsIdx = W == 1;
            yfill = Y;
            yfill(~obsIdx) = nan;
            
            % choose random cluster assignments
            wg.expect.gammank = varargin{1}.wgmm.expect.gammank;
            ridx = argmax(wg.expect.gammank, [], 2);
            
            % go through each cluster and initialize in twostage
            sigmas = zeros(D, D, K);
            mu = zeros(K, D);
            for k = 1:K
                % find mean and fill in means
                mu(k,:) = nanmean(Y(ridx == k, :));
                if ~isclean(mu)
                    sum(ridx == k);
                    warning('mu_k was not clean. Filling in with 0s');
                    mu(k,isnan(mu(k,:))) = 0;
                end

                for i = find(ridx == k)'
                    yfill(i, ~obsIdx(i,:)) = mu(k, ~obsIdx(i,:));
                end
                
                % compute covariance
                sigmas(:,:,k) = cov(yfill(ridx == k, :));

                warning('adding regval in 2S init');
                sigmas(:,:,k) = sigmas(:,:,k) + eye(D) * wg.opts.regularizationValue;
            end
            

            
            % prepare init structure
            wg.params.mu = mu;
            wg.params.sigma = sigmas;
            wg.params.pi = hist(ridx, 1:K) ./ N;
            
        case 'latentMissing-clusterIdx-diagonal'
            % given cluster assignments and initate clusters 
            % in 'diagonal' - Estimate mean, fill NaNs with mean, then estimate diagonal nanvar.
            
            % prepare useful variables
            obsIdx = W == 1;
            yfill = Y;
            yfill(~obsIdx) = nan;
            
            % choose random cluster assignments
            wg.expect.gammank = varargin{1}.wgmm.expect.gammank;
            ridx = argmax(wg.expect.gammank, [], 2);
            
            % go through each cluster and initialize in twostage
            sigmas = zeros(D, D, K);
            mu = zeros(K, D);
            for k = 1:K
                
                % find mean and fill in means
                mu(k,:) = nanmean(Y(ridx == k, :));
                
                % compute covariance
                sigmas(:,:,k) = diag(nanvar(yfill(ridx == k, :)));
            end
            
            % prepare init structure
            wg.params.mu = mu;
            wg.params.sigma = sigmas;
            wg.params.pi = hist(ridx, 1:K) ./ N;
            
        case 'latentSubspace-clusterIdx-twostage'
            % given cluster assignments and initate clusters 
            % in 'twostage' - Estimate mean, fill NaNs with mean, then estimate covar.
            
            % prepare useful variables
            obsIdx = W == 1;
            yfill = Y;
            yfill(~obsIdx) = nan;
            
            % choose random cluster assignments
            wg.expect.gammank = varargin{1}.wgmm.expect.gammank;
            ridx = argmax(wg.expect.gammank, [], 2);
            
            % go through each cluster and initialize in twostage
            sigmas = zeros(D, D, K);
            mu = zeros(K, D);
            for k = 1:K
                % find mean and fill in means
                mu(k,:) = nanmean(Y(ridx == k, :));
                for i = find(ridx == k)'
                    yfill(i, ~obsIdx(i,:)) = mu(k, ~obsIdx(i,:));
                end
                
                % compute covariance
                sigmas(:,:,k) = cov(yfill(ridx == k, :));
                
            end
            
            % prepare init structure
            wg.params.mu = mu;
            wg.params.sigma = sigmas;
            wg.params.pi = hist(ridx, 1:K) ./ N;
            
            % get W, v from sigma.
            for k = 1:K
                [qW(:,:,k), v(k)] = sigma2modelParams(sigmas(:,:,k), wg.opts.model.dopca);
            end
            wg.params.W = qW;
            wg.params.sigmasq = v;
            
        case 'latentSubspace-randW'
            % initialize each cluster by randomly (randn) initializing W, zero means, and random
            % (rand) residual variance
            
            % zero-means
            wg.params.mu = zeros(K, D);
            wg.params.pi = ones(1,K) ./ K;
            wg.params.W = randn(D, wg.opts.model.dopca, K);
            wg.params.sigmasq = rand(1, K);
            
            % initiate sigmas
            for k = 1:K
                W = wg.params.W(:,:,k);
                wg.params.sigma(:,:,k) = W * W' + wg.params.sigmasq(k) .* eye(D);
            end
                
        case 'latentSubspace-randW-mu'
            % initialize each cluster by randomly (randn) initializing W, zero means, and random
            % (rand) residual variance
            
            % zero-means
            wg.params.mu = varargin{1}.wgmm.params.mu;
            wg.params.pi = ones(1,K) ./ K;
            wg.params.W = randn(D, wg.opts.model.dopca, K);
            wg.params.sigmasq = rand(1, K);
            
            % initiate sigmas
            for k = 1:K
                W = wg.params.W(:,:,k);
                wg.params.sigma(:,:,k) = W * W' + wg.params.sigmasq(k) .* eye(D);
            end
            
        case 'latentSubspace-clustIdx'
            % given initial clusters, by randomly (randn) initializing W, zero means, and random
            % (rand) residual variance within each cluster
            
            error('This init doesn''t make much sense now. Other than the mean comptuation, everything is basically still random');
            
            wg.expect.gammank = varargin{1}.wgmm.expect.gammank;
            
            ridx = argmax(wg.expect.gammank, [], 2);
            for k = 1:K
                wg.params.mu(k,:) = nanmean(Y(ridx==k, :), 1);
            end
            wg.params.pi = ones(1,K) ./ K;
            wg.params.W = randn(D, wg.opts.model.dopca, K);
            wg.params.sigmasq = rand(1, K);
            for k = 1:K
                W = wg.params.W(:,:,k);
                wg.params.sigma(:,:,k) = W * W' + wg.params.sigmasq(k) .* eye(D);
            end
            
        case 'latentSubspace-clusterIdx-convergence'
            % given initial cluster guesses, run single-cluster latentSubspace (ppca) on each
            % cluster until convergence, and initialize with these clusters
            wg.expect.gammank = varargin{1}.wgmm.expect.gammank;
            
            ridx = argmax(wg.expect.gammank, [], 2);
            for k = 1:K
                x = Y(ridx == k, :);
                w = W(ridx == k, :);
                d = struct('Y', x, 'W', w, 'K', 1);
                wgi = wgmm(varargin{1}.wgmm.opts, varargin{1}.wgmm.params);
                wgi.params.mu = wgi.params.mu(k, :);
                wgi.params.pi = wgi.params.pi(k);
                wgi.params.W = wgi.params.W(:, :, k);
%                 wgk = wgmmfit(d, 'modelName', 'latentSubspace', 'modelArgs', wg.opts.model, ...
%                     'init', 'latentSubspace-randW-mu', 'verbose', 1, 'replicates', 3, 'MaxIter', 10);
                wgk = wgmmfit(d, 'modelName', 'latentSubspace', 'modelArgs', wg.opts.model, ...
                    'init', 'latentSubspace-randW-mu', 'initArgs', struct('wgmm', wgi), ...
                    'verbose', 1, 'replicates', 3, 'MaxIter', 10);
                %[~, ~, ~, mu, ~, rsltStruct] = ppcax(x, wg.opts.model.dopca, 'Options', struct('Display', 'iter', 'MaxIter', 20, 'TolFun', 1e-3, 'TolX', 1e-3));
                wg.params.mu(k,:) = wgk.params.mu;
                wg.params.W(:,:,k) = wgk.params.W;
%                 if isfield(wgk.params, 'sigma')
%                     wg.params.sigma(:,:,k) = wgk.params.sigma;
%                 end
                wg.params.sigmasq(k) = wgk.params.sigmasq;
                wg.params.pi(k) = sum(ridx == k);
            end
            
        case 'latentSubspace-randn'
            
            wg.params = varargin{1}.wgmm.params;
            
            dLow = wg.opts.model.dopca;
            
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            
            % go through each cluster
            for k = 1:K
                % extract data
                cidx = ridx == k;
                w = W(cidx, :);
                X = Y(cidx, :); 
                
                % compute mu
                wg.params.mu(k, :) = wmean(X, w, 1);
                ws = wg.params.W(:,1:dLow,k);
                newW(:,:,k) = ws .* (1+randn(size(ws)));
            end
            wg.params.W = newW;
            
            
        case 'latentSubspace-clusterIdx-convergence-withrand'

            error('unfinished, but should be quick');
            wg.params = varargin{1}.wgmm.params;
            wg.expect = varargin{1}.wgmm.expect;
            
        case 'latentSubspace-randIdx-convergence'
            % given initial cluster guesses, run single-cluster latentSubspace (ppca) on each
            % cluster until convergence, and initialize with these clusters
            
            ridx = randi([1, K], [N, 1]);
            for k = 1:K
                x = Y(ridx == k, :);
                w = W(ridx == k, :);
                d = struct('Y', x, 'W', w, 'K', 1);
                wgk = wgmmfit(d, 'modelName', 'latentSubspace', 'modelArgs', wg.opts.model, ...
                    'init', 'latentSubspace-randW', 'verbose', 1, 'replicates', 1);
                wg.params.mu(k,:) = wgk.params.mu;
                wg.params.W(:,:,k) = wgk.params.W;
                wg.params.sigma(:,:,k) = wgk.params.sigma;
                wg.params.sigmasq(k) = wgk.params.sigmasq;
                wg.params.pi(k) = sum(ridx == k);
            end
            
            
        case 'latentSubspace-model3'
            % given initial clusters, by randomly (randn) initializing W, zero means, and random
            % (rand) residual variance within each cluster
            assert(isclean(Y), 'latentSubspace-model3 needs DS Y input.');
            
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            
            % go through each cluster
            for k = 1:K
                % extract data
                cidx = ridx == k;
                w = W(cidx, :);
                X = Y(cidx, :); 
                
                % compute mu
                wg.params.mu(k, :) = wmean(X, w, 1);
                xCent = bsxfun(@minus, X, wg.params.mu(k, :));
                
                % model 3 sigma
                wxCent = w .* xCent;
                top = wxCent' * wxCent;
                wtsum = w' * w;
                s = top ./ wtsum;
                
                % decide how to deal with low weight areas
                switch varargin{1}.sigmaCorr
                    case 'zero'
                        s(wtsum < varargin{1}.sigmaCorrThr) = 0;
                    case 'mean'
                        m = mean(s(wtsum(:) > varargin{1}.sigmaCorrThr));
                        s(wtsum < varargin{1}.sigmaCorrThr) = m;
                    case 'diag'
                        s = diag(diag(s));
                    case 'ds'
                        q = (wtsum ./ sum(cidx) ./ mean(w(:)));
                        q = min(q, 1);
                        sc = cov(xCent);
                        % s(wtsum < varargin{1}.sigmaCorrThr) = sc(wtsum < varargin{1}.sigmaCorrThr);
                        s(isnan(s)) = 0;
                        s = s .* q + sc .* (1-q);
                    case 'ds-diag'
                        s = cov(xCent) + diag(diag(s));
                    case 'hack'
                        reconMethod = varargin{1}.sigmaCorrRecon; % usually 'greedy1'
                        mergeMethod = varargin{1}.sigmaCorrMerge; % usually 'wfact-mult-adapt'
                        
                        args = varargin{1}.sigmaCorrMergeArgs;
                        
                        sp = s;
                        sp(isnan(s)) = 0;
                        sr = wgmm.sigmarecon(sp, wtsum, reconMethod);
                        s = wgmm.sigmamerge(sp, sr, wtsum, mergeMethod, args{:});
                    otherwise
                        error('unknown sigma correction method');
                end
                assert(isclean(s))
                
                % assign sigma
                [wg.params.W(:,:,k), wg.params.sigmasq(k)] = sigma2modelParams(s, wg.opts.model.dopca);
                wg.params.pi(k) = sum(ridx==k)./numel(ridx);
            end
            
            
        case 'latentSubspace-iterds'
            initArgs = varargin{1};
            dN = numel(initArgs.dopcas);
            assert(dN >= 1);
            dHigh = size(data.Y, 2);
            
            % dopca has to be adjusted in three palces: 
            % in the initial pca arguments (x2) and in *this* wg options.
            wg.opts.model.dopca = initArgs.dopcas(end);
            
            if dN == 1
                % run normal ds.
                initArgs.wgmm.opts.model.dopca = initArgs.dopcas;
                wginit = initArgs.wgmm;
                
                W = zeros(size(Y, 2), wg.opts.model.dopca, K);
                v = zeros(1, K);
                for k = 1:K
                    [W(:,:,k), v(k)] = sigma2modelParams(wginit.params.sigma(:,:,k), wg.opts.model.dopca);
                    wginit.params.sigma(:,:,k) = W(:,:,k) * W(:,:,k)' + v(k) * eye(dHigh);
                end
                wginit.params.W = W;
                wginit.params.sigmasq = v;
                
            else % recursion
                initArgs.dopcas = initArgs.dopcas(1:(end-1));
                
                wginit = wgmmfit(data, ...
                    'modelName', 'latentSubspace', 'modelArgs', wg.opts.model, ...
                    'init', 'latentSubspace-iterds', 'initArgs', initArgs, ...
                    'verbose', wg.opts.verbose, 'replicates', 1, 'MaxIter', wg.opts.maxIter, ...
                    'TolFun', wg.opts.TolFun);
                
                %yrecon = wginit.recon(data, 'latentMissing');
                %save(sprintf('yrecon%d', wg.opts.model.dopca), 'yrecon');
                
                
                % prepare this init.
                switch initArgs.growW
                    case 'recompute-recon'
                        yrecon = wginit.recon(Y, W, 'latentMissing');
                        mi = argmax(wginit.expect.gammank,[], 2);
                        
                        W = zeros(size(Y, 2), wg.opts.model.dopca, K);
                        v = zeros(1, K);
                        for k = 1:K
                            yy = yrecon(mi == k, :);
                            wginit.params.sigma(:, :, k) = cov(yy);
                            [W(:,:,k), v(k)] = sigma2modelParams(cov(yy), wg.opts.model.dopca);
                        end
                        wginit.params.W = W;
                        wginit.params.sigmasq = v;
                    
                    case 'recompute'
                        W = zeros(size(Y, 2), wg.opts.model.dopca, K);
                        v = zeros(1, K);
                        warning('this is the same as adding 0s I thikn');
                        for k = 1:K
                            oW = wginit.params.W(:, :, k);
                            wginit.params.sigma(:, :, k) = oW * oW' + wginit.params.sigmasq(k) * eye(dHigh);
                            [W(:,:,k), v(k)] = sigma2modelParams(wginit.params.sigma(:, :, k), wg.opts.model.dopca);
                        end
                        wginit.params.W = W;
                        wginit.params.sigmasq = v;
                        
                    case 'add0'
                        disp('add0');
                        [dHigh, oldDLow, K] = size(wginit.params.W);
                        wginit.params.W = [wginit.params.W, zeros(dHigh, wg.opts.model.dopca - oldDLow, K)];
                        
                    case 'addI'
                        disp('addI');
                        [dHigh, oldDLow, K] = size(wginit.params.W);
                        ei = repmat(eye(dHigh), [1,1,K]);
                        wginit.params.W = [wginit.params.W, ei(:, oldDLow+1:wg.opts.model.dopca, :)];
                        
                    case 'addRand'
                        disp('addRand');
                        [dHigh, oldDLow, K] = size(wginit.params.W);
                        wginit.params.W = [wginit.params.W, randn(dHigh, wg.opts.model.dopca - oldDLow, K)];
                        assert(isclean(wginit.params.W));
                        
                    case 'nochange' % don't change w. Should work if coded carefully.
                        error('should be same effect as add0');
                    otherwise
                        error('unknown frow method');
                end
            end
            wg.params = wginit.params;
            wg.expect = wginit.expect;
            wg.stats(1).wginit = wginit;
            
            
        otherwise
            error('wgmm.init: unknown init method');
    end
    
    % check cleanliness and sizes.
    assert(isclean(wg.params.mu));
    if isfield(wg.params, 'sigma')
        assert(isclean(wg.params.sigma));
    end
    assert(isclean(wg.params.pi));
    % assert(all(size(wg.params.mu) == [K, D]), 'The size of init mu is incorrect');
    % assert(all(size(permute(wg.params.sigma, [3 1 2])) == [K, D, D]), 'The size of init sigma is incorrect');
    assert(numel(wg.params.pi) == K, 'The size of init pi is incorrect');
end
