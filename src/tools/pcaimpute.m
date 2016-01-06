function Xi = pcaimpute(X, K, varargin)
% jointly impute missing data and compute PCA via iterative procedure: 
%   estimate PCA and project data
%
%   X - the data, nExperiments x nDims
%       if X has nans, initialization will be done based on nanmean by default
%   K - the number of components to project onto in the second step of the inner loop
%
% Param/Values: 
%   valinit - a full X (with no nans), acting as initialization
%   maxIter - maximum iterations to run
%   trueData - ground truth X (for printing useful error data)
%   realErrFun - error function to compute different between estimate and ground truth 
%   method (for 'ppca' vs 'pcaiter')
%   verbose - logical on weather to print updates. 
%
%   TODO: add "pcainit" to initialize coords directly? useful if use ppca.
%
% Algorithm, see: 
% Ilin and Raiko, 2010. 
%   Practical Approaches to Principal Component Analysis in the Presence of Missing Values
%   Section 2.5
%
% draft code.

    % input checking and initialization
    narginchk(2, inf);
    [X, K, missing, opts] = parseInputs(X, K, varargin{:});
    Xi = X;
    
    % display debug/error information
    if opts.verbose, 
        fprintf('%5s\t%10s\t%10s\t%10s\n', ...
            'rep', 'msd(change)', 'errfun', sprintf('%%expl(K=%d)', K)); 
        fprintf('%5d\t%10.3e\t%10.3e\t%10.2f\n', 0, ...
            msd(Xi(missing), Xi(missing)), opts.realErrFun(Xi, opts.trueData), nan);
    end

    % repeat
    for i = 1:opts.maxIter
        % Step 1.  do PCA
        [coeff, scores, ~, ~, expl, mu] = pca(Xi);
        
        % Step 2. Fill in new missing values
        Xi2 = pcax.recon(coeff(:, 1:K), scores(:, 1:K)', mu')';
        Xip = Xi;
        Xi(missing) = Xi2(missing);
        
        % debug/error information
        if opts.verbose, 
            fprintf('%5d\t%10.3e\t%10.3e\t%10.2f\n', i, ...
                msd(Xip(missing), Xi(missing)), opts.realErrFun(Xi, opts.trueData), sum(expl(1:K)));
        end
    end
    
end

function [X, K, missing, opts] = parseInputs(X, K, varargin)


    p = inputParser();
    p.addRequired('X', @ismatrix);
    p.addRequired('K', @isscalar);
    p.addParameter('valinit', [], @(b) all(size(b) == size(X))); 
    p.addParameter('maxIter', 10, @isscalar);
    p.addParameter('method', 'pcaiter', @(x) sum(strcmp(x, {'pcaiter', 'ppca'})) == 1);
    p.addParameter('verbose', 10, @isscalar);
    p.addParameter('trueData', [], @(b) all(size(b) == size(X)));
    p.addParameter('realErrFun', @(x, y) msd(x(:), y(:)), @isfunc);
    p.parse(X, K, varargin{:});
    
    % get missing data
    missing = isnan(X);
    assert(sum(missing(:)) > 0, 'no missing values detected. nothing for method to impute');
    
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
    assert(~strcmp(p.Results.method, 'ppca'), 'ppca not implemented yet in pcaimpute');
end
