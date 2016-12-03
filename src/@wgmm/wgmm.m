classdef wgmm < handle
    % WGMM: Gaussian Mixture Model with weighted or missing data
    %
    % Behaves similarly to a Gaussian Mixture Model (in fact, in some cases it will behave exactly
    % like a GMM), but allows for the specification of "weights" of each feature in each point. If
    % the weights are logical (0/1), this is treated as a missing data problem. This allows the
    % exploration for a richer set of problems, and leads to different updates than those for the
    % usual GMM.
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
        % wgmm options (controlled by the user)
        %   model (struct) - a struct with model-specific parameters
        %   model.name (chat) - the model used.
        %   init (struct) - optional initialization arguments
        %   init.method (char) - the initialization method
        %   regularizationValue (double) - sigma regularizer
        %   maxIter (int) - maximum number of iterations
        %   replicates (int) - number of replicates to run
        %   TolFun (double) - the stopping tolerance
        %   verbose (logical) - whether to run debug methods
        
        %
        %       sigmaopt (cell) - PROBABLY A MODEL-SPECIFIC OPT. ONLY NEC FOR MODEL3?
        %       (regularizationWeight). 
        %       covarMergeMethod ?
        opts;   

        % expectations (e-step updates)
        % generally, cluster membership will be stored here.
        expect; 
        
        % parameters (m-step updates)
        % usually, cluster mean (mu) and covariances (sigma) and prior (pi)
        % but other model-specific parameters might be stored as well.
        params; 
        
        % general statistics
        stats;  
    end
    
    properties (Hidden)
        
        % helpful implementation structures
        mem;
        % logdetow; % log(det(diag((1./W))));   only needs to be computed once.
        % logdetw;  % log(det(diag((W))));      only needs to be computed once.
    end
    
    methods
        % constructor.
        %   wg = wgmm(opts) initiate options
        %   wg = wgmm(opts, params) include initial params
        function wg = wgmm(opts, varargin)
            if nargin >= 1
                wg.opts = opts;
            else
                wg.opts = wgmm.optionDefaults;
            end
            
            if nargin >= 2
                wg.params = varargin{1}; 
                if isfield(wg.params, 'sigma')
                    assert(size(wg.params.mu, 1) == size(wg.params.sigma, 3));
                end
                assert(size(wg.params.mu, 1) == numel(wg.params.pi));
            end
        end
        
        % general function
        print(wg)
        visualize(wg, varargin);
        
        % EM functions
        logpin = logpost(wg, varargin);
        [ll, expect, wg] = estep(wg, X, W);
        ll = logp(wg, varargin);
        params = mstep(wg, X, W, K);
        wg = init(wg, X, W, K, varargin);
        [sampleX, ks] = sample(gmm, N, X, cidx);
        
        wg = recluster(wg);
        
        
        params = wgmm.mstepModel0(wgmm, data);
        params = wgmm.mstepModel1(wgmm, data);
        params = wgmm.mstepModel3(wgmm, data);
        params = wgmm.mstepModel4exp(wgmm, data);
        params = wgmm.mstepModel5(wgmm, data);
        params = wgmm.mstepLatentSubspace(wgmm, data);
        params = wgmm.mstepLatentMissing(wgmm, data);
        params = wgmm.mstepLatentMissingR(wgmm, data);
    end
    
    methods (Static)
        opts = optionDefaults()
        
        % helper functions
        fwgmm = fit(X, W, K, varargin);
        
        % tools related to log multivariate pdf, specific to wgmm since we work with sigma being
        % DxDxK instead of the usual DxD or DxDxN 
        logp = logmvnpdf(x, mu, sigma, sigmainv);
        
        % helper functions for computing and correcting sigma.
        [sigma, sigmainv, sigmacore, sigmarecon, sigmamerge] = sigmafull(mu, X, W, K, mthods, opts, wg);
        sigmar = sigmarecon(sigma, wtw, method);
        sigma = sigmacore(mu, X, W, K, gammank, coremethod, coreargs, wg);
        sigma = sigmamerge(sigmac, sigmar, wtw, method, varargin);
        
        % general helpful functions
        gmm = gmdist2wgmm(gmdist, X, W, K);
        S = posterior2assignment(post, forcefullassignment);
        S = wgmm2Start(gmm, starttype, varargin)
        compareMeans(gmms, patchSize, titles);
        compare(gmms, patchSize, gmmMethods, varargin);
        compareSigmas(gmms, titles); 
        sm = model5exp(s, X, W, muk, dfn)  
    end
end
