function subvol2ppcawgmmWisoasn(dsSubvolMat, isoSubvolMat, wtSubvolMat, ...
    clusterIdxMat, wgmmMat, isogmmMat, iniFilename, instr)
% dsSubvolMat - matfile name of subvolume
% wtSubvolMat - matfile name of weights
% iniFilename - ini filename with parameters (input)
% clusterIdxMat - matfile name of cluster assignments (if don't have, where to save)
% wgmmMat - matfile name of output wgmm (output)

    params = ini2struct(iniFilename); 

    % rename the parameters
    gmmK = params.gmmK; 
    superPatchSize = params.superPatchSize; 
    patchSize = params.patchSize; 
    ds = params.ds; 
    smallUs = params.smallUs; 
    subtractMean = params.subtractMean;
    minNonNan = params.minNonNan; 
    nPatchesThresh = params.nPatchesThresh; 
    maxECMIter = params.maxECMIter; 
    tolerance = params.tolerance; 
    threshold = params.threshold;
    gmmTolerance = params.gmmTolerance; 
    RECOMPUTECLUSTERS = params.RECOMPUTECLUSTERS; 
    regstring = params.regstring;
    
    nDims = numel(patchSize);

    if exist('instr', 'var')
        fprintf('evaling %s', instr);
        eval(instr);
    end

    %% load subvolumes
    tic
    q = load(isoSubvolMat, 'subvolume'); 
    isoSubvols = q.subvolume;
    q = load(dsSubvolMat, 'subvolume'); 
    dsSubvols = q.subvolume;
    q = load(wtSubvolMat, 'subvolume'); 
    wtSubvols = q.subvolume;
    clear q;
    fprintf('took %5.3f to load the subvolumes\n', toc);
    
    volSize = size(dsSubvols);
    nSubj = volSize(nDims + 1);
    diffPad = (superPatchSize-patchSize)/2;
    assert(any(isIntegerValue(diffPad)), 'difference between superPatchSize and patchSize should be even');
    
    
    %% blur subvolume and train       
    selsub = cropVolume(isoSubvols, [diffPad+1 1], [size(isoSubvols(:,:,:,1))-diffPad nSubj]);
    isoPatches = patchlib.vol2lib(selsub, [patchSize 1]);
     gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', gmmTolerance);

    isoGmmClust = gmdistribution.fit(isoPatches, gmmK, regstring, 1e-4, 'replicates', 1, 'Options', gmmopt);

    postVal = isoGmmClust.posterior(isoPatches);
    [~, clusterIdx] = max(postVal, [], 2);

    % compute the proportion of patches in each cluster
    pis = hist(clusterIdx, 1:gmmK)./numel(clusterIdx); 

    save(clusterIdxMat, 'clusterIdx', 'postVal');
    wg = wgmm(isoGmmClust.mu, isoGmmClust.Sigma, pis);
    save(isogmmMat, 'wg', 'clusterIdx', 'postVal', 'params' ); 

    % compute the proportion of patches in each cluster
    pis = hist(clusterIdx, 1:gmmK)./numel(clusterIdx); 
    
    %% crop patches to 9x9x9
    dsPatches = cropVolume(dsSubvols, [diffPad+1 1], [size(dsSubvols(:,:,:,1))-diffPad nSubj]);
    dsPatches = patchlib.vol2lib(dsPatches, [patchSize 1]);
    
    wtPatches = cropVolume(wtSubvols, [diffPad+1 1], [size(wtSubvols(:,:,:,1))-diffPad nSubj]);
    wtPatches = patchlib.vol2lib(wtPatches, [patchSize 1]);
    
    
    %% sample and initlize covariance/mean

    means = zeros(gmmK, prod(patchSize)); 
    sigmas = zeros(prod(patchSize), prod(patchSize), gmmK);
    
    % sample first.
    nPatches = size(dsPatches,1); 
    randomIdx = randperm(nPatches);
    clusterIdxPerm = clusterIdx(randomIdx); 
    allSamp = []; 
    for k = 1:gmmK
        samp = randomIdx(clusterIdxPerm==k); 
        locsamp = samp(1:min(nPatchesThresh,length(samp)));
        allSamp = [allSamp locsamp]; 
    end
    dsPatches = dsPatches(allSamp,:); 
    wtPatches = wtPatches(allSamp,:); 
    postVal = postVal(allSamp,:); 
    clusterIdx = clusterIdx(allSamp); 
    
    % ppca
    for k = 1:gmmK
        X0 = dsPatches(clusterIdx==k,:);
        assert(~subtractMean, 'subtractMean not well tested yet');
        if subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        W0 = wtPatches(clusterIdx==k,:);       
        
        pk = params.ppcaK;
        d = size(X0, 2);
        tic;
        % compare init costs.
        if params.ecminit
            warning('detected ecminit');
            X0(W0<threshold) = nan; 
            [~, c]  = ecmninitx(X0, 'twostage');
            [u, s, v] = svd(c); 
        else
            [u, s, v] = svd(cov(X0)); 
            X0(W0<threshold) = nan; 
        end
        
        vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
        Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
        opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'iter', 'MaxIter', params.maxPPCAiter);
        [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit);
        if subtractMean
            srecon = bsxfun(@minus, S.Recon, mean(S.Recon, 2));
        else
            srecon = S.Recon;
        end
        means(k,:) = mean(srecon);
        sigmas(:,:,k) = cov(srecon);
        fprintf('ppca cluster %d done in %5.3fs (%d patches)\n', k, toc, size(X0, 1));
        
%         if ~subtractMean
%             % if subtractMean is off, then we did not subtract the 
%             % patchwise mean to create these covariances. However, the
%             % reconstruction code right now assumes we do, and expects mean
%             % and cov as such
%             warning('hacky subtractMean fix. This is very very ugly :(');
%             means(k,:) = mean(bsxfun(@minus, srecon, mean(srecon, 2))); 
%             sigmas(:,:,k) = cov(bsxfun(@minus, srecon, mean(srecon, 2)));
%         end
    end

    %% save wgmm
    wg = wgmm(means, sigmas, pis);
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params' ); 

    warning('remember to add small diagonal component to sigmas in next code'); 


