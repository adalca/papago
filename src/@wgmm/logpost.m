function logpin = logpost(wgmm, data)
% compute log posterior
%   p(params|x) = p(params) * p(x|params) = pi * N(X; mu, owowt .* sigma)
%   logpin (log pi n) is N x K
% method is: 'normal', 'preinvsigma', 'singlesigma' (default)

    % prepare convenient variables
    Y = data.Y;
    if isfield(data, 'W'), wts = data.W; end
    N = size(Y, 1);
    K = size(wgmm.params.mu, 1);
    
    switch wgmm.opts.model.name
        case 'model0' % no weights
             
            logpin = bsxfun(@plus, log(wgmm.params.pi), wgmm.logmvnpdf(Y, wgmm.params.mu, wgmm.params.sigma));
            assert(isclean(logpin), 'log(pi*n) is unclean');
        
        case 'model1-memsafe'
    
            % main implementation, mainly because it's fast.
            % change what sigma is, but need to add back logdetIow.
            logpin = zeros(N, K);
            for k = 1:K

                % first, build owowt
                if isfield(isempty(wgmm.mem, 'logdetow')) && isempty(wgmm.mem.logdetow)
                    % simple and fast. 
                    % log(det(1/D)) = log(\prod(1/D_ii)) = sum(log(1/D_ii)) = sum(-log(D_ii))
                    wgmm.mem.logdetow = sum(-log(wts), 2); 

                    % alternative computation, but this involves a lot of multiplications
                    if wgmm.debug
                        wgmm.mem.logdetow = zeros(N, 1);
                        for i = 1:N
                            wgmm.mem.logdetow(i) = wgmm.logdet(diag(1./wts(i, :)));
                        end
                    end
                    assert(isclean(wgmm.mem.logdetow), 'logdetow is unclean');
                end

                sigma = wgmm.params.sigma(:,:,k);
                Xw = wts .* Y;
                mu = bsxfun(@times, wts, wgmm.params.mu(k, :));
                logpin(:, k) = log(wgmm.params.pi(k)) + wgmm.logmvnpdf(Xw, mu, sigma) - wgmm.mem.logdetow;
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');

            % compute logpost in other methods and test
            if wgmm.debug

                % passing sigma and it needs inverting on the fly. This is very slow.       
                logpinm1 = zeros(N, K);
                for k = 1:K
                    sigma = wgmm.iwAiw(1./wgmm.W, wgmm.params.sigma(:,:,k));
                    logpinm1(:, k) = log(wgmm.params.pi(k)) + wgmm.logmvnpdf(wgmm.X, wgmm.params.mu(k, :), sigma);
                end
                fprintf(2, 'debug: first test method max diff %f\n', max(abs(logpin(:) - logpinm1(:))));

                % pre-compute inversed sigma, which is easy since we have sigmainv and only need to
                % multiply in the W vectors appropriately.
                logpinm2 = zeros(N, K);
                for k = 1:K
                    sigma = wgmm.iwAiw(1./wgmm.W, wgmm.params.sigma(:,:,k));
                    sigmainv = wgmm.iwAiw(wgmm.W, wgmm.params.sigmainv(:,:,k));
                    logpinm2(:, k) = log(wgmm.params.pi(k)) + wgmm.logmvnpdf(wgmm.X, wgmm.params.mu(k, :), sigma, sigmainv);
                end 
                fprintf(2, 'debug: second test method max diff %f\n', max(abs(logpin(:) - logpinm2(:))));
            end
        
        case 'model3'
            logpin = zeros(N, K);            
            
            % first, build logdet(diag(W(i, :)))
            if isempty(wgmm.logdetw)
                % simple and fast.
                % log(det(D)) = log(\prod(D_ii)) = sum(log(D_ii)) = sum(log(D_ii))
                wgmm.logdetw = sum(log(wts), 2);
                
                % alternative computation, but this involves a lot of multiplications
                if wgmm.debug
                    wgmm.logdetw = zeros(N, 1);
                    for i = 1:N
                        wgmm.logdetw(i) = wgmm.logdet(diag(wts(i, :)));
                    end
                end
                assert(isclean(wgmm.logdetw), 'logdetw is unclean');
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
                mu = bsxfun(@times, wts, wgmm.params.mu(k, :)) ;
                if nargin == 1
                    logdetw = wgmm.logdetw;
                    
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

                if ~isempty(wgmm.params.sigmainv)
                    logpin(:, k) = log(wgmm.params.pi(k)) + logmvnpdf(Xw, mu, wgmm.params.sigma(:,:,k), wgmm.params.sigmainv(:,:,k)) - logdetw(:);
                else
                    logpin(:, k) = log(wgmm.params.pi(k)) + logmvnpdf(Xw, mu, wgmm.params.sigma(:,:,k)) - logdetw(:);
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
                muk = wgmm.params.mu(k, :);
                sigmak = wgmm.params.sigma(:,:,k);
                
                for i = 1:N
                    w = wts(i, :);
                    x = Y(i, :);
                    
                    % sigmas
                    Di = wgmm.model4fn(w);
                    sigma = sigmak + Di;
                    
                    % consider pre-computing sigma inverse. But that's slow an dimprecise. 
                    % sigmainv = inv(sigma);
                    logpin(i, k) = log(wgmm.params.pi(k)) + wgmm.logmvnpdf(x, muk, sigma);                    
                end
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');
            
        case {'latentMissing', 'latentSubspace'}
            % missing variables. 
            assert(islogical(wts) | all(wts(:) == 0 | wts(:) == 1));
            
            logmvn = zeros(N, K);
            for i = 1:N
                % extract the observed entry indices for this datapoint
                obsIdx = wts(i, :) == 1;

                % extract the observed data, mu and sigma entries
                yobs = Y(i, obsIdx);
                muobs = wgmm.params.mu(:, obsIdx);
                sigmaobs = wgmm.params.sigma(obsIdx, obsIdx, :);
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                logmvn(i, :) = wgmm.logmvnpdf(yobs, muobs, sigmaobs);
            end
            
            % finally compute the posterior
            logpi = log(wgmm.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
        case {'latentMissingR'}
            
            % extract useful data
            R = data.R;
            ydsmasks = data.ydsmasks;
            yorigs = data.Y;
            
            sigmac = dimsplit(3, wgmm.params.sigma); % split into cell due to faster matlab access
            logmvn = zeros(N, K); tic;
            for i = 1:N
                if mod((i-1), 5000) == 0, fprintf('logpost: %d/%d %3.2fs\n', i, N, toc); tic; end
                
                % extract the observed entry indices for this datapoint
                obsIdx = ydsmasks{i};
                yobs = yorigs{i}(obsIdx);
                r = R.data(R.idx{i}(obsIdx), :); 

                % extract the observed data, mu and sigma entries
                muobs = wgmm.params.mu * r';
                sigmaobs = zeros(sum(obsIdx), sum(obsIdx), K);
                for k = 1:K
                    sigmaobs(:,:,k) = r * sigmac{k} * r';
                end
                
                % compute compute the multivariate normal for each k via logN(y^Oi; mu^Oi, sigma^Oi)
                logmvn(i, :) = wgmm.logmvnpdf(yobs, muobs, sigmaobs);
            end
            fprintf('logpost: %d/%d %3.2fs\n', i, N, toc);
            
            % finally compute the posterior
            logpi = log(wgmm.params.pi);
            logpin = bsxfun(@plus, logpi, logmvn);
            
        otherwise
            error('unknown logp method');
    end
    
    assert(isclean(logpin), 'logp(params|x) is not clean');
end
