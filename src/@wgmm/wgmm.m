classdef wgmm < handle
    % WGMM Weighted GMM class
    %   Weighted Gaussian Mixture Model (WGMM). Behaves similarly to a Gaussian Mixture Model (in
    %   fact, in some cases it will behave exactly like a GMM), but allows for the specification of
    %   "weights" of each feature in each point. This allows the exploration for a richer set of
    %   problems, and leads to slightly different updates than those for the usual GMM.
    %
    % Weighted GMM might go by other names as well, such as sparse GMM.
    %
    % This class offers object features, including fitting of the model, and offers static functions
    % for testing and analysis of these models.
    %
    % Since we allow each feature of each observation to have a weight, we need to deal with the
    % case where some feature in the computation of a cluster mean has no overall weight (TODO), and
    % the case wehre some *pair* of features has no overall weight (i.e. it could be that each
    % of the pair of features has support in general, but they never both have support in the same
    % observation), affecting the cluster covariation (sigma) estimates. We deal with the latter by
    % assuming there is structure in the covariance that allows us to reconstruct the missing
    % pairing.
    
    properties
        % main properties
        mu;
        sigma;
        pi;
        
        % derivative properties, useful to have around, especially during fitting. 
        % This should be optional.
        sigmainv;
    end
    
    properties (Hidden)
        
        % fitting parameters
        initmethod = 'exemplar';
        debug = false;
        sigmareg = nan;
        sigmaopt = {0};
        
        % Update parameters. 
        % Most of these available methods should get cleaned up after development.
        covarUpdateMethod = 'model3';
        muUpdateMethod = 'model3';
        logpUpdateMethod = 'model3';
        
        covarReconMethod = 'greedy1';
        covarMergeMethod = 'wfact-mult-adapt'; % 'none'
        
        % helpful implementation structures
        logdetow; % log(det(diag((1./W)))); only needs to be computed once.
        logdetw; % log(det(diag((W)))); only needs to be computed once.
    end
    
    methods
        % constructor
        function gmm = wgmm(varargin)
            if nargin >= 1, gmm.mu = varargin{1}; end
            if nargin >= 2, gmm.sigma = varargin{2}; end
            if nargin >= 3, gmm.pi = varargin{3}; end
            if nargin >= 4, gmm.sigmainv = varargin{4}; end
                
            assert(size(gmm.mu, 1) == size(gmm.sigma, 3));
            assert(size(gmm.mu, 1) == numel(gmm.pi));
        end
        
        function print(wg)
            s = '--- WGMM object:\n';
            
            s = sprintf('%soperating methods:', s);
            s = sprintf('%s%30s:%s\n', s, 'initialization', wg.initmethod);
            s = sprintf('%s%30s:%s\n', s, 'covariance update', wg.covarUpdateMethod);
            s = sprintf('%s%30s:%s\n', s, 'mu update', wg.muUpdateMethod);
            s = sprintf('%s%30s:%s\n', s, 'logp update', wg.logpUpdateMethod);
            s = sprintf('%s%30s:%s\n', s, 'covariance reconstruction', wg.covarReconMethod);
            s = sprintf('%s%30s:%s\n', s, 'covariance merge', wg.covarMergeMethod);
            fprintf(s);
            sys.structmemory(wg)
        end
        
        % EM functions
        logpin = logpost(wg, varargin);
        [gammank, ll] = estep(wg, X, W);
        ll = logp(wg, varargin);
        wg = mstep(wg, X, W, K, gammank);
        wg = init(wg, X, W, K, varargin);
        visualize(wg, varargin);
        [sampleX, ks] = sample(gmm, N, X, cidx);
    end
    
    methods (Static)
        %% helper functions
        fwgmm = fit(X, W, K, varargin);
        %s = iwAiw(W, A); % TODO: I'm not sure this is used anymore
        %sxr = sx(X, S, dodebug); % TODO: I'm not sure this is used anymore
        
        % tools related to log multivariate pdf, specific to wgmm since we work with sigma being
        % DxDxK instead of the usual DxD or DxDxN 
        logp = logmvnpdf(x, mu, sigma, sigmainv);
        
        % helper functions for computing and correcting sigma.
        [sigma, sigmainv, sigmacore, sigmarecon, sigmamerge] = sigmafull(mu, X, W, K, gammank, mthods, opts);
        sigmar = sigmarecon(sigma, wtw, method);
        sigma = sigmacore(mu, X, W, K, gammank, coremethod, coreargs);
        sigma = sigmamerge(sigmac, sigmar, wtw, method, varargin);
        
        % general helpful functions
        gmm = gmdist2wgmm(gmdist, X, W, K);
        S = posterior2assignment(post, forcefullassignment);
        S = wgmm2Start(gmm, starttype, varargin)
        compareMeans(gmms, patchSize, titles);
        compare(gmms, patchSize, gmmMethods, varargin);
        compareSigmas(gmms, titles); 
        
    end
    
end
