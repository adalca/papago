function [logpin, varargout] = logpost(gmm, X, W)
% compute log posterior
%   p(params|x) = p(params) * p(x|params) = pi * N(X; mu, owowt .* sigma)
%   logpin (log pi n) is N x K
% method is: 'normal', 'preinvsigma', 'singlesigma' (default)

    % prepare convenient variables
    [N, D] = size(X);    
    K = size(gmm.mu, 1);
    varargout = {};
    
    switch gmm.logpUpdateMethod
        case 'model0' % no weights
             
            logpin = bsxfun(@plus, log(gmm.pi), wgmm.logmvnpdf(X, gmm.mu, gmm.sigma));
            assert(isclean(logpin), 'log(pi*n) is unclean');
        
        case 'model1-memsafe'
    
            % main implementation, mainly because it's fast.
            % change what sigma is, but need to add back logdetIow.
            logpin = zeros(N, K);
            for k = 1:K

                % first, build owowt
                if isempty(gmm.logdetow)
                    % simple and fast. 
                    % log(det(1/D)) = log(\prod(1/D_ii)) = sum(log(1/D_ii)) = sum(-log(D_ii))
                    gmm.logdetow = sum(-log(W), 2); 

                    % alternative computation, but this involves a lot of multiplications
                    if gmm.debug
                        gmm.logdetow = zeros(N, 1);
                        for i = 1:N
                            gmm.logdetow(i) = wgmm.logdet(diag(1./W(i, :)));
                        end
                    end
                    assert(isclean(gmm.logdetow), 'logdetow is unclean');
                end

                sigma = gmm.sigma(:,:,k);
                Xw = W .* X;
                mu = bsxfun(@times, W, gmm.mu(k, :));
                logpin(:, k) = log(gmm.pi(k)) + wgmm.logmvnpdf(Xw, mu, sigma) - gmm.logdetow;
            end
            assert(isclean(logpin), 'log(pi*n) is unclean');

            % compute logpost in other methods and test
            if gmm.debug

                % passing sigma and it needs inverting on the fly. This is very slow.       
                logpinm1 = zeros(N, K);
                for k = 1:K
                    sigma = wgmm.iwAiw(1./gmm.W, gmm.sigma(:,:,k));
                    logpinm1(:, k) = log(gmm.pi(k)) + wgmm.logmvnpdf(gmm.X, gmm.mu(k, :), sigma);
                end
                fprintf(2, 'debug: first test method max diff %f\n', max(abs(logpin(:) - logpinm1(:))));

                % pre-compute inversed sigma, which is easy since we have sigmainv and only need to
                % multiply in the W vectors appropriately.
                logpinm2 = zeros(N, K);
                for k = 1:K
                    sigma = wgmm.iwAiw(1./gmm.W, gmm.sigma(:,:,k));
                    sigmainv = wgmm.iwAiw(gmm.W, gmm.sigmainv(:,:,k));
                    logpinm2(:, k) = log(gmm.pi(k)) + wgmm.logmvnpdf(gmm.X, gmm.mu(k, :), sigma, sigmainv);
                end 
                fprintf(2, 'debug: second test method max diff %f\n', max(abs(logpin(:) - logpinm2(:))));
            end
        
        case 'model3'
            logpin = zeros(N, K);            
            
            % first, build logdet(diag(W(i, :)))
            if isempty(gmm.logdetw)
                % simple and fast.
                % log(det(D)) = log(\prod(D_ii)) = sum(log(D_ii)) = sum(log(D_ii))
                gmm.logdetw = sum(log(W), 2);
                
                % alternative computation, but this involves a lot of multiplications
                if gmm.debug
                    gmm.logdetw = zeros(N, 1);
                    for i = 1:N
                        gmm.logdetw(i) = gmm.logdet(diag(W(i, :)));
                    end
                end
                assert(isclean(gmm.logdetw), 'logdetw is unclean');
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
                Xw = X .* W;
                mu = bsxfun(@times, W, gmm.mu(k, :)) ;
                if nargin == 1
                    logdetw = gmm.logdetw;
                    
                else
                    if all(W(:) == 1)
                        logdetw = zeros(size(W, 1), 1);
                    else
                        % compute logdet(diag(W(i, :)); for each i
                        % logdetw = zeros(N, 1);
                        % for i = 1:N
                        %    logdetw(i) = logdet(diag(W(i, :)));
                        % end
                        logdetw = sum(log(W), 2);
                    end
                end
                % Note: passing in the inverse sigma is much faster, but might lose some accuracy.
                % from our tests, at least on usRate = 2, patchSize of 5^3, the maximum error is ~1e-10.

                if ~isempty(gmm.sigmainv)
                    logpin(:, k) = log(gmm.pi(k)) + logmvnpdf(Xw, mu, gmm.sigma(:,:,k), gmm.sigmainv(:,:,k)) - logdetw(:);
                else
                    logpin(:, k) = log(gmm.pi(k)) + logmvnpdf(Xw, mu, gmm.sigma(:,:,k)) - logdetw(:);
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
            murnk = zeros(N, K, D);
            for k = 1:K
                muk = gmm.mu(k, :);
                sigmak = gmm.sigma(:,:,k);
                
                for i = 1:N
                    w = W(i, :);
                    x = X(i, :);
                    
                    % sigmas
                    Di = gmm.model4fn(w);
                    sigma = sigmak + Di;
                    
                    % consider pre-computing sigma inverse. But that's slow an dimprecise. 
                    % sigmainv = inv(sigma);
                    logpin(i, k) = log(gmm.pi(k)) + wgmm.logmvnpdf(x, muk, sigma);                    
                    murnk(i, k, :) = reshape((sigmak / sigma) * (x - muk)' + muk', [1, 1, D]);
                end
            end
            varargout{1} = murnk;
            assert(isclean(logpin), 'log(pi*n) is unclean');
            assert(isclean(murnk), 'log(pi*n) is unclean');
            
        otherwise
            error('unknown logp method');
    end
end
