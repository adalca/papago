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
            wg.expect = varargin{1}.wgmm.expect;
            
            if strcmp(wg.opts.init.method, 'wgmmr')
                wg.params.W = wg.params.W + randn(size(wg.params.W)) * 0.05;
                wg.params.sigma = wg.wv2sigma();
            end
            
            if strcmp(wg.opts.init.method, 'wgmmI')
                e = repmat(eye(D), [1,1,K]);
                wg.params.W = wg.params.W + e(:, 1:size(wg.params.W, 2), :);
                wg.params.sigma = wg.wv2sigma();
            end            
            
            if any(strcmp(wg.opts.model.name, {'latentSubspace', 'wLatentSubspace'})) && ...
                (~isfield(wg.params, 'W') || size(wg.params.W, 2) ~= wg.opts.model.dopca)
                warning('wgmm init: W not the right size. Re-estimating from Sigma');

                W = zeros(size(Y, 2), wg.opts.model.dopca, K);
                v = zeros(1, K);
                for k = 1:K
                    [W(:,:,k), v(k)] = sigma2modelParams(wg.params.sigma(:,:,k), wg.opts.model.dopca);
                    wg.params.sigma(:,:,k) = W(:,:,k) * W(:,:,k)' + v(k) * eye(size(Y, 2));
                    assert(rank(W(:,:,k)) == wg.opts.model.dopca);
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
            
        case 'latentSubspace-wgmm-convergence'
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
                wgi.params.W = wgi.params.W(:, :, k); % + randn(size(wgi.params.W(:, :, k))) * 0.01;
                wgi.params.sigmasq(k) = wgi.params.sigmasq;
                
                wgk = wgmmfit(d, 'modelName', wg.opts.model.name, 'modelArgs', wg.opts.model, ...
                    'init', 'wgmm', 'initArgs', struct('wgmm', wgi), ...
                    'verbose', 1, 'replicates', 1, 'MaxIter', 21);
                
                wg.params.mu(k,:) = wgk.params.mu;
                wg.params.W(:,:,k) = wgk.params.W;
                wg.params.sigmasq(k) = wgk.params.sigmasq;
                wg.params.pi(k) = sum(ridx == k)./numel(ridx);
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
                newW(:,:,k) = ws .* (1+randn(size(ws))*0.2);
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
            
        case 'LS-blurdct'
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            dLow = wg.opts.model.dopca;
            
            blurSigma = varargin{1}.blurSigma;
            patchSize = varargin{1}.patchSize;
            
            stdObs = std(data.Y(data.W(:) == 1));

            
            % go through each cluster
            for k = 1:K
                % get dct 
                nW = dct(eye(D));
                nWc = cellfunc(@(x) reshape(x, patchSize), dimsplit(2, nW));
                % blur the dct components
                nWcb = cellfunc(@(x) volblur(x, blurSigma), nWc);
                nWcs = cellfunc(@(x) x(:) ./ max(abs(x(:))) .* stdObs, nWcb);
                nW = cat(2, nWcs{:});
                
                cidx = ridx == k;
                wg.params.mu(k, :) = wmean(Y(cidx, :), data.W(cidx, :), 1);
                wg.params.W(:,:,k) = nW(:, 1:dLow);
                wg.params.sigmasq(k) = varargin{1}.wgmm.params.sigmasq(k);
                wg.params.pi(k) = mean(ridx == k);
            end
            
            
        case 'LS-perlin'
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            dLow = wg.opts.model.dopca;
            
            blurSigma = varargin{1}.blurSigma;
            patchSize = varargin{1}.patchSize;
            
            stdObs = std(data.Y(data.W(:) == 1));
            
            % go through each cluster
            for k = 1:K
                bs = linspace(blurSigma, 0.5, dLow);
                nWcs = arrayfunc(@(x) volblur(perlin(patchSize), x), bs);
                nWcs = cellfunc(@(x) x - mean(x(:)), nWcs);
                nWcs = cellfunc(@(x) x(:)./ max(abs(x(:))) .* stdObs, nWcs);
                nW = cat(2, nWcs{:});
                
                cidx = ridx == k;
                wg.params.mu(k, :) = wmean(Y(cidx, :), data.W(cidx, :), 1);
                wg.params.W(:,:,k) = nW(:, 1:dLow);
                wg.params.sigmasq(k) = varargin{1}.wgmm.params.sigmasq(k);
                wg.params.pi(k) = mean(ridx == k);
            end
            
        case 'LS-diag-perlin'
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            dLow = wg.opts.model.dopca;
            
            blurSigma = varargin{1}.blurSigma;
            patchSize = varargin{1}.patchSize;
            
            stdObs = std(data.Y(data.W(:) == 1));
            
            % go through each cluster
            for k = 1:K
                ei = eye(D);
                bs = linspace(blurSigma, 0.5, dLow);
                nWcs = arrayfunc(@(x) volblur(perlin(patchSize), x), bs);
                nWcs = cellfunc(@(x) x - mean(x(:)), nWcs);
                nWcs = cellfunc(@(x) x(:)./ max(abs(x(:))) .* stdObs, nWcs);
                nW = cat(2, nWcs{:}) + ei(:, 1:dLow);
                
                cidx = ridx == k;
                wg.params.mu(k, :) = wmean(Y(cidx, :), data.W(cidx, :), 1);
                wg.params.W(:,:,k) = nW(:, 1:dLow);
                wg.params.sigmasq(k) = varargin{1}.wgmm.params.sigmasq(k);
                wg.params.pi(k) = mean(ridx == k);
            end
            
        case 'LS-dct'
            % compute expectations and cluster assignments
            [~, wg.expect] = varargin{1}.wgmm.estep(data); 
            ridx = argmax(wg.expect.gammank, [], 2);
            dLow = wg.opts.model.dopca;
            
            % go through each cluster
            for k = 1:K
                % get dct 
                nW = dct(eye(D));
                cidx = ridx == k;
                wg.params.mu(k, :) = wmean(Y(cidx, :), data.W(cidx, :), 1);
                wg.params.W(:,:,k) = nW(:, 1:dLow);
                wg.params.sigmasq(k) = varargin{1}.wgmm.params.sigmasq(k);
                wg.params.pi(k) = mean(ridx == k);
            end
            
            
        case 'latentSubspace-model3'
            % given initial clusters, by randomly (randn) initializing W, zero means, and random
            % (rand) residual variance within each cluster
            assert(isclean(Y), 'latentSubspace-model3 needs DS Y input.');
            
            % compute expectations and cluster assignments
            if isfield(data, 'allW')
                d = struct('Y', data.allY, 'W', data.allW, 'K', data.K);
                [~, wg.expect] = varargin{1}.wgmm.estep(d); 
                ridx = argmax(wg.expect.gammank, [], 2);                
            else
                [~, wg.expect] = varargin{1}.wgmm.estep(data); 
                ridx = argmax(wg.expect.gammank, [], 2);
            end
            
            % go through each cluster
            for k = 1:K
                % extract data
                cidx = ridx == k;
                if isfield(data, 'allW')
                    w = data.allW(cidx, :);
                    X = data.allY(cidx, :);
                elseif isfield(data, 'lowW')
                    w = data.lowW(cidx, :);
                    X = Y(cidx, :); 
                else
                    w = data.W(cidx, :);
                    X = Y(cidx, :); 
                end
                
                
                % compute mu
                wg.params.mu(k, :) = wmean(X, w, 1);
                xCent = bsxfun(@minus, X, wg.params.mu(k, :));
                
                % model 3 sigma
                wxCent = w .* xCent;
                top = wxCent' * wxCent;
                wtsum = w' * w;
                s = top ./ wtsum;
                
%                 warning('testing with iso');
%                 s = cov(data.isodata(cidx, :));
%                 wtsum = wtsum*0+1000;
                
%                 warning('testing with ds');
%                 s = cov(xCent);
%                 wtsum = wtsum*0+1000;
                
                % decide how to deal with low weight areas
                switch varargin{1}.sigmaCorr
                    case 'zero'
                        s(wtsum < varargin{1}.sigmaCorrThr) = 0;
                    case 'mean'
                        m = mean(s(wtsum(:) > varargin{1}.sigmaCorrThr));
                        s(wtsum < varargin{1}.sigmaCorrThr) = m;
                    case 'diag'
                        s = diag(diag(s));
                    case 'rand'
                        thr = varargin{1}.sigmaCorrThr;
                        badidx = wtsum < thr | isnan(s);
                        r = normrnd(0, std(s(~badidx)), size(s));
                        tmap = triu(true(size(s)));
                        s(badidx & tmap) = r(badidx & tmap);
                        r = r';
                        s(badidx & ~tmap) = r(badidx & ~tmap);
                    case 'flip-rand'
                        thr = varargin{1}.sigmaCorrThr;
                        badidx = wtsum < thr;
                        
                        % rotate along 1st dim
                        xCent1 = rot90_3D_stack(xCent, varargin{1}.patchSize, 1);
                        w1 = rot90_3D_stack(w, varargin{1}.patchSize, 1);
                        s1 = ( (w1 .* xCent1)' * (w1 .* xCent1) ) ./ (w1' * w1);
                        s(badidx) = s1(badidx);
                        badidx = badidx & ( (w1' * w1) < thr);
                        s(badidx) = nan;

                        % rotate along 2st dim
                        xCent1 = rot90_3D_stack(xCent, varargin{1}.patchSize, 2);
                        w1 = rot90_3D_stack(w, varargin{1}.patchSize, 2);
                        s1 = ( (w1 .* xCent1)' * (w1 .* xCent1) ) ./ (w1' * w1);
                        s(badidx) = s1(badidx);
                        badidx = badidx & ( (w1' * w1) < thr);
                        s(badidx) = nan;
                        
                        % rotate along 3st dim
                        xCent1 = rot90_3D_stack(xCent, varargin{1}.patchSize, 3);
                        w1 = rot90_3D_stack(w, varargin{1}.patchSize, 3);
                        s1 = ( (w1 .* xCent1)' * (w1 .* xCent1) ) ./ (w1' * w1);
                        s(badidx) = s1(badidx);
                        badidx = badidx & ( (w1' * w1) < thr);
                        s(badidx) = nan;
                        
                        r = normrnd(0, std(s(~badidx)), size(s));
                        tmap = triu(true(size(s)));
                        s(badidx & tmap) = r(badidx & tmap);
                        r = r';
                        s(badidx & ~tmap) = r(badidx & ~tmap);
                    case 'ds-randn'
                        badidx = wtsum < varargin{1}.sigmaCorrThr;
                        sc = cov(xCent) + randn(size(s)) * std(s(~badidx(:)));
                        s(badidx) = sc(badidx);
                    case 'ds'
%                         q = (wtsum ./ sum(cidx) ./ mean(w(:)));
%                         q = min(q, 1);
                        sc = cov(xCent);
                        s(wtsum < varargin{1}.sigmaCorrThr) = sc(wtsum < varargin{1}.sigmaCorrThr);
%                         s(isnan(s)) = 0;
%                         s = s .* q + sc .* (1-q);
                    case 'ds-diag'
                        s = cov(xCent) + diag(diag(s));
                    case 'hack'
                        reconMethod = varargin{1}.sigmaCorrRecon; % usually 'greedy1'
                        mergeMethod = varargin{1}.sigmaCorrMerge; % usually 'wfact-mult-adapt'
                        
                        args = varargin{1}.sigmaCorrMergeArgs;
                        
                        sp = s;
                        sp(wtsum(:) < varargin{1}.sigmaCorrThr) = 0;
%                         varargin{1}.sigmaCorrThr
                        wtw = wtsum;
                        wtw(wtsum(:) < varargin{1}.sigmaCorrThr) = 0;
                        sr = wgmm.sigmarecon(sp, wtw, reconMethod);
                        s = wgmm.sigmamerge(sp, sr, wtsum, mergeMethod, args{:});
                    otherwise
                        error('unknown sigma correction method');
                end
                assert(isclean(s))
                
                % assign sigma
                [wg.params.W(:,:,k), wg.params.sigmasq(k)] = sigma2modelParams(s, wg.opts.model.dopca);
                wg.params.pi(k) = mean(ridx==k);
               
                if isfield(varargin{1}, 'sigmaCorrConv') && varargin{1}.sigmaCorrConv
                    if isfield(data, 'allW')
                        d = struct('Y', data.allY(cidx, :), 'W', data.allW(cidx, :), 'K', 1);
                    else
                        d = struct('Y', data.Y(cidx, :), 'W', data.W(cidx, :), 'K', 1);
                    end
                    warning('change w to real w not low...?');
                    wgi = wgmm(varargin{1}.wgmm.opts, varargin{1}.wgmm.params);
                    wgi.params.mu = wg.params.mu(k, :);
                    wgi.params.pi = 1;
                    wgi.params.W = wg.params.W(:, :, k); % + randn(size(wgi.params.W(:, :, k))) * 0.01;
                    wgi.params.sigmasq(k) = wg.params.sigmasq(k);

                    wgk = wgmmfit(d, 'modelName', wg.opts.model.name, 'modelArgs', wg.opts.model, ...
                        'init', 'wgmm', 'initArgs', struct('wgmm', wgi), ...
                        'verbose', 1, 'replicates', 1, 'MaxIter', 5);
                    
                    wg.params.mu(k,:) = wgk.params.mu;
                    wg.params.W(:,:,k) = wgk.params.W;
                    wg.params.sigmasq(k) = wgk.params.sigmasq;
                    wg.params.pi(k) = mean(ridx == k);
                end
                
            end
            
            
        case {'latentSubspace-iterds', 'latentSubspace-iterds-zeromean', 'latentSubspace-iterds-randv'}
            initArgs = varargin{1};
            wgz = initArgs.wgmm;
            initMethod = wg.opts.init.method;
            dN = numel(initArgs.dopcas);
            assert(dN >= 1);
            % dHigh = size(data.Y, 2);
            
            wginit_prev = [];
            
            for di = 1:numel(initArgs.dopcas)
                pca = initArgs.dopcas(di);
                
                % init params
                params = wgz.params;
                
                % get parameters for new run from prev params
                params.W = []; params.sigmasq = [];
                for k = 1:K
                    [params.W(:,:,k), params.sigmasq(k)] = sigma2modelParams(params.sigma(:,:,k), pca);
                end
                
                if strcmp('latentSubspace-iterds-zeromean', initMethod)
                    fprintf('zeroing mean ...');
                    mu = mu * 0;
                end

                if strcmp('latentSubspace-iterds-randv', initMethod)
                    fprintf('setting v0 to rand...');
                    params.sigmasq = rand(1, K);
                end
                
                mi = argmax(wgz.expect.gammank, [], 2);
                for k = 1:K
                    kdata = struct('Y', data.Y(mi == k, :), 'W', data.W(mi == k, :), 'K', 1);
                    kparams = struct('mu', params.mu(k, :), 'W', params.W(:,:,k), 'pi', 1, 'sigmasq', params.sigmasq(k));
                    kwginit = wgmm(initArgs.wgmm.opts, kparams);
                    wgzk = wgmmfit(kdata, ...
                        'modelName', 'latentSubspace', 'modelArgs', struct('dopca', pca), ...
                        'init', 'wgmm', 'initArgs', struct('wgmm', kwginit), ...
                        'verbose', wg.opts.verbose, 'replicates', 1, 'MaxIter', 2, ...
                        'MinIter', 1, 'TolFun', wg.opts.TolFun);
                    
                    params.mu(k, :) = wgzk.params.mu;
                    params.W(:, :, k) = wgzk.params.W;
                    params.sigmasq(k) = wgzk.params.sigmasq;
                end
                
                wginit = wgmm(initArgs.wgmm.opts, params);
                wginit.expect = wgz.expect;
                
                wgz = wgmmfit(data, ...
                    'modelName', 'latentSubspace', 'modelArgs', struct('dopca', pca), ...
                    'init', 'wgmm', 'initArgs', struct('wgmm', wginit), ...
                    'verbose', wg.opts.verbose, 'replicates', 1, 'MaxIter', wg.opts.maxIter, ...
                    'MinIter', wg.opts.minIter, 'TolFun', wg.opts.TolFun);
                wgz.stats(1).wginit = wginit_prev;
                wginit_prev = wgz;
                
                assert(~isfield(wgz.params, 'sigma'))
                wgz.params.sigma = wgz.wv2sigma;
            end
            
            % do nothing.
            wg.params = wgz.params;
            wg.expect = wgz.expect;
            wg.stats = wgz.stats;
            wg.opts.maxIter = 0;
            wg.opts.minIter = -1;

%             
%             % dopca has to be adjusted in three palces: 
%             % in the initial pca arguments (x2) and in *this* wg options.
%             wg.opts.model.dopca = initArgs.dopcas(end);
%             
%             if dN == 1
%                 % run normal ds.
%                 initArgs.wgmm.opts.model.dopca = initArgs.dopcas;
%                 wginit = initArgs.wgmm;
%                 
%                 W = zeros(size(Y, 2), wg.opts.model.dopca, K);
%                 v = zeros(1, K);
%                 for k = 1:K
%                     [W(:,:,k), v(k)] = sigma2modelParams(wginit.params.sigma(:,:,k), wg.opts.model.dopca);
%                     wginit.params.sigma(:,:,k) = W(:,:,k) * W(:,:,k)' + v(k) * eye(dHigh);
%                 end
%                 wginit.params.W = W;
%                 wginit.params.sigmasq = v;
%                 
%             else % recursion
%                 initArgs.dopcas = initArgs.dopcas(1:(end-1));
%                 
%                 wginit = wgmmfit(data, ...
%                     'modelName', 'latentSubspace', 'modelArgs', wg.opts.model, ...
%                     'init', wg.opts.init.method, 'initArgs', initArgs, ...
%                     'verbose', wg.opts.verbose, 'replicates', 1, 'MaxIter', wg.opts.maxIter, ...
%                     'MinIter', wg.opts.minIter, 'TolFun', wg.opts.TolFun);
%                 
%                 %yrecon = wginit.recon(data, 'latentMissing');
%                 %save(sprintf('yrecon%d', wg.opts.model.dopca), 'yrecon');
%                 
%                 % prepare this init.
%                 switch initArgs.growW
%                     case 'recompute-recon'
%                         yrecon = wginit.recon(Y, W, 'latentMissing');
%                         mi = argmax(wginit.expect.gammank,[], 2);
%                         
%                         W = zeros(size(Y, 2), wg.opts.model.dopca, K);
%                         v = zeros(1, K);
%                         for k = 1:K
%                             yy = yrecon(mi == k, :);
%                             wginit.params.sigma(:, :, k) = cov(yy);
%                             [W(:,:,k), v(k)] = sigma2modelParams(cov(yy), wg.opts.model.dopca);
%                         end
%                         wginit.params.W = W;
%                         wginit.params.sigmasq = v;
%                     
%                     case 'recompute'
%                         W = zeros(size(Y, 2), wg.opts.model.dopca, K);
%                         v = zeros(1, K);
%                         warning('this is the same as adding 0s I thikn');
%                         for k = 1:K
%                             oW = wginit.params.W(:, :, k);
%                             wginit.params.sigma(:, :, k) = oW * oW' + wginit.params.sigmasq(k) * eye(dHigh);
%                             [W(:,:,k), v(k)] = sigma2modelParams(wginit.params.sigma(:, :, k), wg.opts.model.dopca);
%                         end
%                         wginit.params.W = W;
%                         wginit.params.sigmasq = v;
%                         
%                     case 'add0'
%                         disp('add0');
%                         [dHigh, oldDLow, K] = size(wginit.params.W);
%                         wginit.params.W = [wginit.params.W, zeros(dHigh, wg.opts.model.dopca - oldDLow, K)];
%                         
%                     case 'addI'
%                         disp('addI');
%                         [dHigh, oldDLow, K] = size(wginit.params.W);
%                         ei = repmat(eye(dHigh), [1,1,K]);
%                         wginit.params.W = [wginit.params.W, ei(:, oldDLow+1:wg.opts.model.dopca, :)];
%                         
%                     case 'addRand'
%                         disp('addRand');
%                         [dHigh, oldDLow, K] = size(wginit.params.W);
%                         wginit.params.W = [wginit.params.W, randn(dHigh, wg.opts.model.dopca - oldDLow, K)*0.001];
%                         assert(isclean(wginit.params.W));
%                         
%                     case 'nochange' % don't change w. Should work if coded carefully.
%                         error('should be same effect as add0');
%                     otherwise
%                         error('unknown frow method');
%                 end
%             end
%             wginit.params = rmfield(wginit.params, 'sigma');
%             wg.params = wginit.params;
%             wg.expect = wginit.expect;
%             wg.stats(1).wginit = wginit;
            
            
            
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


function Y = rot90_3D(X, norotdim)

    Y = cellfunc(@squeeze, dimsplit(norotdim, X));
    Y = cellfunc(@rot90, Y);

    % put back together
    switch norotdim
        case 1
            Y = cellfunc(@(v) permute(v, [3, 1, 2]), Y);
        case 2
            Y = cellfunc(@(v) permute(v, [1, 3, 2]), Y);
        case 3
            Y = cellfunc(@(v) permute(v, [1, 2, 3]), Y);
    end
       
    Y = cat(norotdim, Y{:});
end

function Z = rot90_3D_stack(xCent, patchSize, norotdim)
    Xc = dimsplit(1, xCent);
    Xc = cellfunc(@(x) reshape(x, patchSize), Xc);

    % rotate along 1st dim
    Z = cellfunc(@(x) rot90_3D(x, norotdim), Xc);
    Z = cellfunc(@(x) x(:)', Z);
    Z = cat(1, Z{:});
end
