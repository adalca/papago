function logpin = logpost(wg, data)
% compute log posterior
%   p(params|x) = p(params) * p(x|params) = pi * N(X; mu, owowt .* sigma)
%   logpin (log pi n) is N x K
% method is: 'normal', 'preinvsigma', 'singlesigma' (default)

    % prepare convenient variables
    Y = data.Y;
    if isfield(data, 'W'), wts = data.W; end
    if isfield(data, 'estepW') && numel(wg.stats) < 10
        warning('estep hack (top)'); 
        wts = data.estepW == 1; 
    end
    [N, dHigh] = size(Y);
    K = size(wg.params.mu, 1);
    
    switch wg.opts.model.name
        case 'model0' % no weights
             
            logpin = bsxfun(@plus, log(wg.params.pi), wgmm.logmvnpdf(Y, wg.params.mu, wg.params.sigma));
            assert(isclean(logpin), 'log(pi*n) is unclean');
        
        case 'model1-memsafe'
    
            % main implementation, mainly because it's fast.
            % change what sigma is, but need to add back logdetIow.
            logpin = zeros(N, K);
            for k = 1:K

                % first, build owowt
                if isfield(isempty(wg.mem, 'logdetow')) && isempty(wg.mem.logdetow)
                    % simple and fast. 
                    % log(det(1/D)) = log(\prod(1/D_ii)) = sum(log(1/D_ii)) = sum(-log(D_ii))
                    wg.mem.logdetow = sum(-log(wts), 2); 

                    % alternative computation, but this involves a lot of multiplications
                    if wg.debug
                        wg.mem.logdetow = zeros(N, 1);
                        for i = 1:N
                            wg.mem.logdetow(i) = wg.logdet(diag(1./wts(i, :)));
                        end
                    end
                    assert(isclean(wg.mem.logdetow), 'logdetow is unclean');
                end

                sigma = wg.params.sigma(:,:,k);
                Xw = wts .* Y;
                mu = bsxfun(@times, wts, wg.params.mu(k, :));
                logpin(:, k) = log(wg.params.pi(k)) + wgmm.logmvnpdf(Xw, mu, sigma) - wg.mem.logdetow;
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');

            % compute logpost in other methods and test
            if wg.debug

                % passing sigma and it needs inverting on the fly. This is very slow.       
                logpinm1 = zeros(N, K);
                for k = 1:K
                    sigma = wg.iwAiw(1./wg.W, wg.params.sigma(:,:,k));
                    logpinm1(:, k) = log(wg.params.pi(k)) + wgmm.logmvnpdf(wg.X, wg.params.mu(k, :), sigma);
                end
                fprintf(2, 'debug: first test method max diff %f\n', max(abs(logpin(:) - logpinm1(:))));

                % pre-compute inversed sigma, which is easy since we have sigmainv and only need to
                % multiply in the W vectors appropriately.
                logpinm2 = zeros(N, K);
                for k = 1:K
                    sigma = wg.iwAiw(1./wg.W, wg.params.sigma(:,:,k));
                    sigmainv = wg.iwAiw(wg.W, wg.params.sigmainv(:,:,k));
                    logpinm2(:, k) = log(wg.params.pi(k)) + wgmm.logmvnpdf(wg.X, wg.params.mu(k, :), sigma, sigmainv);
                end 
                fprintf(2, 'debug: second test method max diff %f\n', max(abs(logpin(:) - logpinm2(:))));
            end
        
        case 'model3'
            logpin = zeros(N, K);            
            
            % first, build logdet(diag(W(i, :)))
            if isempty(wg.logdetw)
                % simple and fast.
                % log(det(D)) = log(\prod(D_ii)) = sum(log(D_ii)) = sum(log(D_ii))
                wg.logdetw = sum(log(wts), 2);
                
                % alternative computation, but this involves a lot of multiplications
                if wg.debug
                    wg.logdetw = zeros(N, 1);
                    for i = 1:N
                        wg.logdetw(i) = wg.logdet(diag(wts(i, :)));
                    end
                end
                assert(isclean(wg.logdetw), 'logdetw is unclean');
            end
            
            for k = 1:K

                % how to actually compute? 
                % method 1: just do as if no weights involved
                % method 2: somehow "count" the high weights more in the probability. 
                %   this would maybe be achieved by collapsing those dimensions? so maybe
                %   (x-mu) --> w .* (x - mu); and
                %   sigma --> w .* w' * sigma. 
                %   but this computation leads to cancelling of ws in the logmvnpdf. 
                %   Sounds Counterintuitive :(
                %   just need to offset by logdetw.
                % see 09/14/2015 onenote discussion.
                sys.warn(['model 3 logpost method not fully decided\n', ...
                    'Option 1: W on (X-mu) only, not sigma: low-W dims matter less. But leads to bad assignments in many cases(?)', ...
                    'Option 2: just normal X, mu, sigma. Perhaps problem if initial X is really bad?', ...
                    'Currently: Option 1'], ...
                    'SingleWarn', true);
                Xw = Y .* wts;
                mu = bsxfun(@times, wts, wg.params.mu(k, :)) ;
                if nargin == 1
                    logdetw = wg.logdetw;
                    
                else
                    if all(wts(:) == 1)
                        logdetw = zeros(size(wts, 1), 1);
                    else
                        % compute logdet(diag(W(i, :)); for each i
                        % logdetw = zeros(N, 1);
                        % for i = 1:N
                        %    logdetw(i) = logdet(diag(W(i, :)));
                        % end
                        logdetw = sum(log(wts), 2);
                    end
                end
                % Note: passing in the inverse sigma is much faster, but might lose some accuracy.
                % from our tests, at least on usRate = 2, patchSize of 5^3, the maximum error is ~1e-10.

                if ~isempty(wg.params.sigmainv)
                    logpin(:, k) = log(wg.params.pi(k)) + logmvnpdf(Xw, mu, wg.params.sigma(:,:,k), wg.params.sigmainv(:,:,k)) - logdetw(:);
                else
                    logpin(:, k) = log(wg.params.pi(k)) + logmvnpdf(Xw, mu, wg.params.sigma(:,:,k)) - logdetw(:);
                end
                    
%                 sigma = gmm.sigma(:,:,k);
%                 Xw = W .* X;
%                 mu = bsxfun(@times, W, gmm.mu(k, :));
%                 for i = 1:size(X, 1)
%                     w = W(i, :)';
%                     s = (w * w') .* sigma;
%                     s = sigma;
%                     x = Xw(i, :);
%                     m = mu(i, :);                    
%                     try
%                         logpin(i, k) = log(gmm.pi(k)) + gmm.logmvnpdf(x, m, s);
%                     catch
%                         logpin(i, k) = log(gmm.pi(k)) + gmm.logmvnpdf(x, m, s + eye(size(s,1)) * 0.1);
%                     end
%                 end
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');
            
        case {'model4', 'model5'}
            % compute gammank. This will be slow since we have to invert for each * subject. 
            % perhaps we could approximate it?
            % also compute mu^r_nk
            logpin = zeros(N, K);
            for k = 1:K
                muk = wg.params.mu(k, :);
                sigmak = wg.params.sigma(:,:,k);
                
                for i = 1:N
                    w = wts(i, :);
                    x = Y(i, :);
                    
                    % sigmas
                    Di = wg.model4fn(w);
                    sigma = sigmak + Di;
                    
                    % consider pre-computing sigma inverse. But that's slow an dimprecise. 
                    % sigmainv = inv(sigma);
                    logpin(i, k) = log(wg.params.pi(k)) + wgmm.logmvnpdf(x, muk, sigma);                    
                end
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');
            
        case 'latentMissing'
            % missing variables. 
            assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
            
            logmvn = zeros(N, K);
            for i = 1:N
                % extract the observed entry indices for this datapoint
                obsIdx = wts(i, :) == 1;

                % extract the observed data, mu and sigma entries
                yobs = Y(i, obsIdx);
                muobs = wg.params.mu(:, obsIdx);
                sigmaobs = wg.params.sigma(obsIdx, obsIdx, :);
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                % lmvn = wgmm.logmvnpdf(yobs, muobs, sigmaobs); % older and slower
                lmvn = experimentalLogmvnpdf(yobs, muobs, sigmaobs); % new method, using ECMOBJ basically.
                logmvn(i, :) = lmvn;
            end
            
            % finally compute the posterior
            logpi = log(wg.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
        case {'latentSubspace', 'gpuLatentSubspace'}
            
            if strcmp(wg.opts.model.name, 'gpuLatentSubspace')
                Y = data.gpuY;
                wts = data.gpuW;
                wg.params.mu = gpuArray(single(wg.params.mu));
                wg.params.W = gpuArray(single(wg.params.W));
                
                warning('hack: e-step is not gpu-ized');
                Y = gather(Y);
                wts = gather(wts);
                wg.params.mu = gather(wg.params.mu);
                wg.params.W = gather(wg.params.W);
                
                logmvn = zeros(N, K, class(Y));
            else
                % missing variables. 
                assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1), 'Weights are non binary');
                logmvn = zeros(N, K);
            end
            
            
           
            
            for i = 1:N
                % extract the observed entry indices for this datapoint
                obsIdx = wts(i, :) == 1;

                % extract the observed data, mu and sigma entries
                yobs = Y(i, obsIdx);
                muobs = wg.params.mu(:, obsIdx);
                % sigmaobs = wg.params.sigma(obsIdx, obsIdx, :);
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                % lmvn = wgmm.logmvnpdf(yobs, muobs, sigmaobs); % older and slower
                %lmvn = experimentalLogmvnpdf(yobs, muobs, sigmaobs); % new method, using ECMOBJ basically.
                Wobs = wg.params.W(obsIdx, :, :);
                lmvn = experimentalLogmvnpdfW(yobs, muobs, Wobs, wg.params.sigmasq);
                logmvn(i, :) = lmvn;
            end
            
            % finally compute the posterior
            logpi = log(wg.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
        case {'wLatentSubspace'}
            % missing variables. 
            
            logmvn = zeros(N, K);
            for i = 1:N
                % extract the observed entry indices for this datapoint
                obsIdx = wts(i, :) > 0.5;

                % extract the observed data, mu and sigma entries
                yobs = Y(i, obsIdx);
                muobs = wg.params.mu(:, obsIdx);
                % sigmaobs = sigma(obsIdx, obsIdx, :);
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                % lmvn = wgmm.logmvnpdf(yobs, muobs, sigmaobs); % older and slower
                %lmvn = experimentalLogmvnpdf(yobs, muobs, sigmaobs); % new method, using ECMOBJ basically.
                Wobs = wg.params.W(obsIdx, :, :);
                lmvn = experimentalLogmvnpdfW(yobs, muobs, Wobs, wg.params.sigmasq);
                logmvn(i, :) = lmvn;
            end
            
            % finally compute the posterior
            logpi = log(wg.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
        case 'latentSubspaceR'
            warning('Still unsure of how to treat R and G at the edges. Need to fix this!');
            
            % extract useful data
            R = data.R;
            ydsmasks = data.ydsmasksFullVoxels; % only use the "absolutely" full voxels.
            yorigs = data.Y;
            
            % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows
            % -- so we transpose the R data matrix and index in columns instead of rows.
            RdataTrans = R.data';

            % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than
            % indexing into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times
            % per iteration, it ends up mattering for our iterations.
            wCell = dimsplit(3, wg.params.W); % split into cell due to faster matlab access
            logmvn = zeros(N, K); tic;
            for i = 1:N
                % adding a bit of verbosity
                if wg.opts.verbose >= 2 && mod((i-1), 5000) == 0, 
                    fprintf('logpost: %d/%d %3.2fs\n', i, N, toc); tic; 
                end
                
                % extract the observed entry indices for this datapoint
                obsIdx = ydsmasks{i};
                yobs = yorigs{i}(obsIdx);
                r = RdataTrans(:, R.idx{i}(obsIdx))';
                
                % extract the observed data, mu and sigma entries
                muobs = wg.params.mu * r';
                
                wobs = cell(K, 1);
                for k = 1:K
                    wobs{k} = r * wCell{k};
                end
                
                lmvn = experimentalLogmvnpdfW(yobs, muobs, wobs, wg.params.sigmasq); 
                logmvn(i, :) = lmvn;
            end
            if wg.opts.verbose >= 2
                fprintf('logpost: %d/%d %3.2fs\n', i, N, toc);
            end
            
            % finally compute the posterior
            logpi = log(wg.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
                
        case 'latentMissingR'
            warning('Still unsure of how to treat R and G at the edges. Need to fix this!');
            warning('LMR E-Step: doing a trick or rotating data back in atlas space :(');
            
            % extract useful data
            R = data.R;
            ydsmasks = data.ydsmasksFullVoxels; % only use the "absolutely" full voxels.
            yorigs = data.Y;
            
            % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows
            % -- so we transpose the R data matrix and index in columns instead of rows.
            RdataTrans = R.data';

            % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than
            % indexing into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times
            % per iteration, it ends up mattering for our iterations.
            sigmaCell = dimsplit(3, wg.params.sigma); % split into cell due to faster matlab access
            logmvn = zeros(N, K); tic;
            for i = 1:N
                % adding a bit of verbosity
                if wg.opts.verbose >= 2 && mod((i-1), 5000) == 0, 
                    fprintf('logpost: %d/%d %3.2fs\n', i, N, toc); tic; 
                end
                
                % extract the observed entry indices for this datapoint
                obsIdx = ydsmasks{i};
                yobs = yorigs{i}(obsIdx);
                r = RdataTrans(:, R.idx{i}(obsIdx))';
                
                % extract the observed data, mu and sigma entries
                muobs = wg.params.mu * r';
                sigmaobs = cell(K, 1);
                for k = 1:K
                    sigmaobs{k} = r * sigmaCell{k} * r';
                end
                
                % hack
%                 g = data.G.data(:, data.G.idx{i});
%                 ww = data.ydsmasks{i} * g';
%                 obsIdx = ww > 0.5;
%                 yobs = yorigs{i} * g(obsIdx, :)';
%                 muobs = wg.params.mu(:, obsIdx);
%                 sigmaobs = wg.params.sigma(obsIdx, obsIdx, :);
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                % lmvn = wgmm.logmvnpdf(yobs, muobs, sigmaobs); % older and slower
                lmvn = experimentalLogmvnpdf(yobs, muobs, sigmaobs); % new method, using ECMOBJ basically.
                logmvn(i, :) = lmvn;
            end
            if wg.opts.verbose >= 2
                fprintf('logpost: %d/%d %3.2fs\n', i, N, toc);
            end
            
            % finally compute the posterior
            logpi = log(wg.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
%             figure(); imagesc(logpin)
            
        otherwise
            error('unknown logp method');
    end
    
    assert(isclean(logpin), 'logp(params|x) is not clean');
end

function obj = experimentalLogmvnpdf(yobs, muobs, sigmaobs)
% using the tricks from the ecmobj function. should only be used when we
% know there's missing values, otherwise it's unnecessarily slow

    isc = iscell(sigmaobs);

    if strcmp(class(yobs), 'gpuArray')
        obj = zeros(1, size(muobs, 1), 'single', 'gpuArray');
    else
        obj = zeros(1, size(muobs, 1));
    end
    for k = 1:size(muobs, 1)

        ltp = log(2*pi);
        obj(k) = -0.5 * numel(yobs) * ltp;

        if isc
            [SubChol, CholState] = chol(sigmaobs{k});
        else
            [SubChol, CholState] = chol(sigmaobs(:,:,k));
        end

        if CholState > 0
            error(message('finance:ecmnobj:NonPosDefSubCovar'));
        end

        SubResid = SubChol' \ (yobs(:) - muobs(k, :)');

        obj(k) = obj(k) - 0.5 * (SubResid' * SubResid);
        obj(k) = obj(k) - sum(log(diag(SubChol)));
    end
end

function obj = experimentalLogmvnpdfW(yobs, muobs, Wobs, v)

    isc = iscell(Wobs);
    if isc
        K = numel(Wobs);
        [dHigh, dLow] = size(Wobs{1});
    else
        [dHigh, dLow, K] = size(Wobs);
    end
    ltp = log(2*pi);

    if isa(yobs, 'gpuArray')
        obj = zeros(1, size(muobs, 1), 'single', 'gpuArray') - 0.5 * numel(yobs) * ltp;
    else
        obj = zeros(1, size(muobs, 1)) - 0.5 * numel(yobs) * ltp;
    end
    
    % for each cluster;
    for k = 1:K
        % extrtact W and L
        if isc, W = Wobs{k}; else, W = Wobs(:,:,k); end
        L = W ./ sqrt(v(k));
        df = (yobs(:) - muobs(k, :)');
        sLow = L' * L + eye(dLow);
        
        % method 1
        if dLow < 30 % 30 is estimated...
            q = df' * L;
            t0 = q * (sLow \ q');
            t1 = 0.5 * (df' * df  - t0) ./ v(k);
            % via Sylvester's determinant identity (SDI)
            % det(s) = det(W * W' + v I
            %        = det(v * L * L' + v I
            %        = det(v (L * L' + I)
            %        = v^dHigh * det(L' * L + I) % By SDI
            t2 = 0.5 * (dHigh * log(v(k))) + 0.5 * logdet(sLow);
        else
        
            % method 2
            SubChol = chol(sLow);
            t0p1 = SubChol' \ L';
            t0 = t0p1' * (t0p1 * df); % do the df mult first!
            t1 = 0.5 * df' * (df - t0) ./ v(k);
            t2 = 0.5 * (dHigh * log(v(k))) + sum(log(diag(SubChol)));
        end
        
        obj(k) = obj(k) - t1 - t2;
    end
end

