function [gmm, patches] = learngmm(volumes, params, opt, wopt)
% LEARNGMM learn gmms from (sub) volumes
%
%   [gmm, patches] = learngmm(volumes, K, popt, wopt) train weighted gaussian mixture models based
%   on the given volumes. Volumes is a nVolumes-by-1 or nVolumes-by-1
%
%   params:
%       patchSize
%       K
%       nSamples
%       minWeight
%
%   popt: (papago options)
%       subtractpatchmean, 
%       globalinit, <globalgmm>, 
%       debugisogmm
%       
%   wopt: (wgmm options)
%       updateMethod (model - we likely want model0 (unweighted) or model3)
%       replicates
%       TolFun
%       regularizationValue
%       regularizationWeight

    % TODO: work with subvol2XW

    % should only have X or X,W
    assert(ismember(size(volumes, 2), [1, 2]));
    K = params.K;

    % prepare sampling amounts. 
    if params.nSamples == 0 % all
        nSamples = numel(volumes) * numel(volumes{1});
    elseif params.nSamples < 1 % fraction of all
        nSamples = numel(volumes) * numel(volumes{1}) * params.nSamples;
    else 
        nSamples = params.nSamples;
    end
    
    % prepare volumes
    if size(volumes, 2) > 1
        v = dimsplit(2, volumes);
    else
        v = {volumes};
    end
    
    % sample patches
    patches = patchlib.vol2samples(nSamples, params.patchSize, v{:});
    
    % prepare data (X,W)
    if size(volumes, 2) == 2
        X = patches{1};
        W = max(patches{2}, params.minWeight);
    else
        X = patches;
        W = ones(size(X));
    end

    % subtract patch means if necessary
    if opt.subtractpatchmean;
        X = bsxfun(@minus, X, mean(X, 2));
    end
    
    % prepare global-based initialization, if necessary
    if opt.globalinit
        % re-evaluate weights
        pgmm = opt.globalgmm;

        % get the posteriors of the local patches. 
        % pgmm is expected to be a matlab gmdist for now. use posterior()
        post = pgmm.posterior(X); 
        spost = sum(post); % TODO: make sure this is not pre-normalized? or it's ok if it is?

        % keep the highest K clusters
        [~, si] = sort(spost, 'descend'); 
        clusts = si(1:K);
        spost = spost(clusts);

        pi = spost ./ sum(spost);
        gmm = wgmm(X, X*0+1, K, pgmm.mu(clusts, :), pgmm.Sigma(:,:,clusts), ones(1, K)./K);
        q = gmm.estep();
        pi = sum(q) ./ sum(q(:));
        gmm = wgmm(X, X*0+1, K, pgmm.mu(clusts, :), pgmm.Sigma(:,:,clusts), pi);
        wopt.init = gmm;
    end
    
    % run the new GMM. 
    gmm = wgmmfit(X, W, K, wopt);
    
    % compute a matlab-gmm as well (this should only be done on iso data)
    if opt.debugisogmm 
        gmdist = fitgmdist(X, K, opt.debugisogmmargs{:});
    end
