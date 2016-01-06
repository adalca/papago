function [Xi, gmd] = gmmimpute(X, nClust, K, varargin)
% jointly impute missing data and compute PCA.
%
% Param/Values: valinit, maxIter, trueData, realErrFun, method (for ppca vs pcaiter), verbose, 
%   TODO: also "pcainit" to initialize coords directly? useful if use ppca.
%
% draft code.
%
% Major TODO, implement tipping & bigshop mixture of (p)pca 1997/99
%
% X is of size nExperiments x nDims

    % input checking and initialization
    narginchk(2, inf);
    [X, nClust, K, missing, opts] = parseInputs(X, nClust, K, varargin{:});
    Xi = X;
    
    % display debug/error information
    if opts.verbose
        fprintf('%5s\t%10s\t%10s\t%10s\n', 'rep', 'msd(change)', 'errfun', 'negloglik'); 
        fprintf('%5d\t%3.2e\t%3.2e\t%3.2f\n', 0, ...
            msd(Xi(missing), Xi(missing)), opts.realErrFun(Xi, opts.trueData), nan);
    end

    % repeat
    for i = 1:opts.maxIter
        % Step 1.  do GMM
        st = statset('Display', 'off', 'MaxIter', 5);
        gmd = fitgmdist(Xi, nClust, 'Replicates', 5, 'Regularize', 0.00001, 'Options', st);
        
        if ~isempty(opts.valwts)
            [~, ci] = max(posterior(gmd, Xi, opts.valwts), [], 2); % get cluster assignment
        else
            [~, ci] = max(posterior(gmd, Xi), [], 2);
        end
        
        % Step 2. Fill in new missing values using inpaint
        XiDebug = Xi; assert(isclean(XiDebug));
        
        for c = 1:nClust
            cidx = find(ci == c);
            
            k = min(K, sum(cidx));
            if k == 0
                warning('Found a cluster with no points :(');
                continue;
            end
            if k < K, warning('cluster %d has only %d points\n', c, k); end
            
            % TODO: imputation methods should be:
            % - just inpaint with covariances
            % - pca
            % - pcaimpute
            % - ppca
            
            % impute without pca. This is using out own implementation of covariance-based
            % inpainting
            if K == 0
                Xm = X; Xm(missing) = nan;
                Xt = gmmInpaint(Xm, gmd, ci, c);
                Xi(cidx, :) = Xt(cidx, :);
                err('todo: instead of hard assignment in missing, soft assignment with weights!');
                % see implementation below.
                
            else
                tmpXi = Xi(cidx, :);
            
                
                % impute with pca
                [coeff, scores, latent, t2, expl, mu] = pca(tmpXi);
                if k > size(coeff, 2)
                    warning('coeffs is %dx%d, k=%d (K=%d)', size(coeff), k, K);
                    k = size(coeff, 2);
                end
                if k == 0
                    warning('Found a cluster with no points :(');
                    continue;
                end
                Xirecon = pcax.recon(coeff(:, 1:k), scores(:, 1:k)', mu)';
                warning('ok, but this is basically using the original values :( the gmm doesn''t do that much');
                
                % re-assign values
                if ~isempty(opts.valwts) % if weights are given
                    Xi(cidx, :) = X(cidx, :) .* opts.valwts(cidx, :) + Xirecon .* (1 - opts.valwts(cidx, :));
                else
                    tmpXi(missing(cidx, :)) = Xirecon(missing(cidx, :));
                    Xi(cidx, :) = tmpXi;
                end
                
               
                
                % impute with pcaimpute. % on a quick try, this actually seems worse. prob since
                % only initializing within the cluster (?)
%                 a = X(cidx, :); m = missing(cidx, :); a(m) = nan;
%                 piopts = {'trueData', opts.trueData(cidx, :), 'realErrFun', opts.realErrFun, 'maxIter', 10};
%                 Xi(cidx, :) = pcaimpute(a, K, 'valinit', tmpXi, piopts{:});
%                 opts.realErrFun(Xi, opts.trueData)

    
                % try ppca
                % st = statset('ppca'); st.Display = 'final'; st.MaxIter = 250;
                % [coeff,score,pcvar, mu, v, S] = ppca(X(cidx, :), k, 'Options', st);
                % fprintf('ppca err: %f\n', opts.realErrFun(S.Recon, opts.trueData(cidx, :)))
            end
        end
        assert(isclean(Xi));
        
        if opts.verbose
            fprintf('%5d\t%3.2e\t%3.2e\t%3.2f', i, ...
                msd(XiDebug(missing), Xi(missing)), opts.realErrFun(Xi, opts.trueData), gmd.NegativeLogLikelihood);
        end
        fprintf('%3.2f ', gmd.ComponentProportion);
        fprintf('\tHist:');
        h = hist(ci, 1:nClust);
        fprintf('%3.2f ', h/sum(h));
        fprintf('\n');
    end
end


function [X, nClust, K, missing, opts] = parseInputs(X, nClust, K, varargin)

    p = inputParser();
    p.addRequired('X', @ismatrix);
    p.addRequired('nClust', @isscalar);
    p.addRequired('K', @isscalar);
    p.addParameter('valinit', [], @(b) all(size(b) == size(X))); 
    p.addParameter('valwts', [], @(b) all(size(b) == size(X))); 
    p.addParameter('maxIter', 10, @isscalar);
    p.addParameter('pcamethod', 'pcaiter', @(x) sum(strcmp(x, {'pcaiter', 'ppca'})) == 1);
    p.addParameter('verbose', 10, @isscalar);
    p.addParameter('trueData', [], @(b) all(size(b) == size(X)));
    p.addParameter('realErrFun', @(x, y) msd(x(:), y(:)), @isfunc);
    p.parse(X, nClust, K, varargin{:});
    
    % get missing data
    missing = isnan(X);
    if isempty(p.Results.valwts)
        assert(sum(missing(:)) > 0, 'no missing values detected. nothing for method to impute');
    end
    
    % initialize X
    if ~isempty(p.Results.valinit)
        cond = X(~missing) == p.Results.valinit(~missing);
        assert(all(cond(:)));
        X = p.Results.valinit;
    else
        meanvalrep = repmat(nanmean(X), [size(X, 1), 1]);
        X(missing) = meanvalrep(missing);
    end
    assert(isclean(X));

    % options (mostly for maxIter, method, verbose, trueData, realErrFun
    opts = p.Results;
    
    % ppca not implemented yet :(
    assert(~strcmp(p.Results.pcamethod, 'ppca'), 'ppca not implemented yet in pcaimpute');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMM Functions to be able to compute weighted posterior computations...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [post, NlogL] = posterior(gmd, X, W)

    if nargin <= 2
        [post, NlogL] = gmd.posterior(X);
    else
        % this is the standard gmd.posterior, but calling wdensity with a weight. see "avd"
        % modification in our wdensity() function below.
        covNames = { 'diagonal','full'};
        CovType = find(strncmpi(gmd.CovType,covNames,length(gmd.CovType)));
        log_lh = wdensity(X, gmd.mu, gmd.Sigma, gmd.PComponents, gmd.SharedCov, CovType, W);
        [ll, post] = gmdestep(log_lh);
        NlogL = -ll;
    end
end

function [ll, post] = gmdestep(log_lh)
    maxll = max (log_lh,[],2);
    %minus maxll to avoid underflow
    post = exp(bsxfun(@minus, log_lh, maxll));
    %density(i) is \sum_j \alpha_j P(x_i| \theta_j)/ exp(maxll(i))
    density = sum(post,2);
    %normalize posteriors
    post = bsxfun(@rdivide, post, density);
    logpdf = log(density) + maxll;
    ll = sum(logpdf) ;
end

function   [log_lh,mahalaD]=wdensity(X, mu, Sigma, p, sharedCov, CovType, W)
%WDENSITY Weighted conditional density and mahalanobis distance.
%   LOG_LH = WDENSITY(...) returns log of component conditional density
%   (weighted by the component probability) of X. LOG_LH is a N-by-K matrix
%   LOG_LH, where K is the number of Gaussian components. LOG_LH(I,J) is
%   log (Pr(point I|component J) * Prob( component J))
%
%   [LOG_LH, MAHALAD]=WDENSITY(...) returns the Mahalanobis distance in
%   the N-by-K matrix MAHALAD. MAHALAD(I,J) is the Mahalanobis distance of
%   point I from the mean of component J.

%   Copyright 2007 The MathWorks, Inc.


    log_prior = log(p);
    [n,d]=size(X);
    k=size(mu,1);
    log_lh = zeros(n,k);
    mahalaD = zeros(n,k);
    logDetSigma = -Inf;
    for j = 1:k
        if sharedCov
            if j == 1
                if CovType == 2 % full covariance
                    [L,f] = chol(Sigma);
                    diagL = diag(L);
                    if (f ~= 0)|| any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                        error(message('stats:gmdistribution:wdensity:IllCondCov'));
                    end
                    logDetSigma = 2*sum(log(diagL));
                else %diagonal
                    L = sqrt(Sigma);
                    if  any(L < eps( max(L))*d)
                          error(message('stats:gmdistribution:wdensity:IllCondCov'));
                    end
                    logDetSigma = sum( log(Sigma) );
                end
            end
        else %different covariance
            if CovType == 2 %full covariacne
                % compute the log determinant of covariance
                [L,f] = chol(Sigma(:,:,j) );
                diagL = diag(L);
                if (f ~= 0) || any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                     error(message('stats:gmdistribution:wdensity:IllCondCov'));
                end
                logDetSigma = 2*sum(log(diagL));
            else %diagonal covariance
                L = sqrt(Sigma(:,:,j)); % a vector
                if  any(L < eps(max(L))*d)
                     error(message('stats:gmdistribution:wdensity:IllCondCov'));
                end
                logDetSigma = sum( log(Sigma(:,:,j)) );

            end
        end

        % start avd modification
        if exist('W', 'var')
            Xcentered = W .* bsxfun(@minus, X, mu(j,:));
        else
            Xcentered = bsxfun(@minus, X, mu(j,:));
        end
        % end avd modification
       
        if CovType == 2
            xRinv = Xcentered /L ;
        else
            xRinv = bsxfun(@times,Xcentered , (1./ L));
        end
        mahalaD(:,j) = sum(xRinv.^2, 2);

        log_lh(:,j) = -0.5 * mahalaD(:,j) +...
            (-0.5 *logDetSigma + log_prior(j)) - d*log(2*pi)/2;
        %get the loglikelihood for each point with each component
        %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
    end
end
   

