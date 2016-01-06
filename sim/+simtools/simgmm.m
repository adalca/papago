function [X, mu, sigma, pi, cidx, varargout] = simgmm(N, D, K, method, varargin)
% simulate data for gmm
% methods: 'simple' or 'frompatchesgmm'


    varargout = {};

    switch method
        case 'simple'
            if numel(varargin) > 0
                musepfact = varargin{1}; % how separate the clusters should be
            else
                musepfact = 10; % how separate the clusters should be
            end

            % simulate random generative parameters
            mu = bsxfun(@times, rand(K, D), musepfact*(1:K)');
            sigma = zeros(D, D, K);
            sigmainv = zeros(D, D, K);
            for i = 1:K
                a = randn(D);
                sigma(:,:,i) = a' * a;
                sigmainv(:,:,i) = inv(sigma(:,:,i));
            end
            pi = rand(1, K); pi = pi ./ sum(pi);
            pi = pi * 0 + 1/K; % equal

            % create data: sample, no extra noise
            cidx = discretesample(pi, N);
            X = zeros(N, D);
            for i = 1:K
                X(cidx == i, :) = mvnrnd(mu(i, :)', sigma(:, :, i), sum(cidx == i));
            end
            
            varargout{1} = sigmainv;
            
        case 'frompatchesgmm'
            patches = varargin{1};
            
            % run gmm on given patches. patches should be [NxD]
            assert(all(size(patches, 2) == D));
            
            st = statset('Display','iter');
            if ispc
                gmm = fitgmdist(patches, K, 'RegularizationValue', 0.00000001, 'replicates', 3, 'Options', st);
                pi = gmm.ComponentProportion;
            else
                warning('upgrade to MATLAB 2015 and change this back ASAP...Katie'); 
                gmm = fitgmdist(patches, K, 'Regularize', 0.00000001, 'replicates', 3);
                pi = gmm.PComponents; 
            end
            
            mu = gmm.mu;
            sigma = gmm.Sigma;
            
            % create data: sample, no extra noise
            cidx = discretesample(pi, N);
            X = zeros(N, D);
            for i = 1:K
                X(cidx == i, :) = mvnrnd(mu(i, :)', sigma(:, :, i), sum(cidx == i));
            end
            
            sigmainv = zeros(D, D, K);
            for i = 1:K
                sigmainv(:,:,i) = inv(sigma(:,:,i));
            end
            
            varargout{1} = sigmainv;
            
        otherwise
            error('unknown sumlation method');
            
    end
