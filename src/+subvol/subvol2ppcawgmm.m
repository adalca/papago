function subvol2ppcawgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
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
    threshold = params.threshold; 
    subtractMean = params.subtractMean;
    minNonNan = params.minNonNan; 
    nPatchesThresh = params.nPatchesThresh; 
    maxECMIter = params.maxECMIter; 
    tolerance = params.tolerance; 
    gmmTolerance = params.gmmTolerance; 
    RECOMPUTECLUSTERS = params.RECOMPUTECLUSTERS; 
    regstring = params.regstring;
    
    nDims = numel(patchSize);

    %% load subvolumes
    tic
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
    if (~exist(clusterIdxMat, 'file') || RECOMPUTECLUSTERS)
        
        sigmaBlur = 1/3 * (ds / smallUs);
   
        % weight-based blurring
        dsBlurvol = volblur(dsSubvols .* wtSubvols, sigmaBlur); 
        wtBlurvol = volblur(wtSubvols, sigmaBlur); 
        smallBlurPatches = dsBlurvol./wtBlurvol; 
        dsBlurvol = []; % clear volumes
        wtBlurvol = []; % clear volumes
        
        % get library of blurred patches.
        smallBlurPatches = cropVolume(smallBlurPatches, [diffPad+1 1], [size(dsSubvols(:,:,:,1))-diffPad nSubj]);
        smallBlurPatches = patchlib.vol2lib(smallBlurPatches, [patchSize 1]);
        
        % resize to small blurred patches
        nPatches = size(smallBlurPatches,1); 
        subsampledSize = (ceil((ceil(smallUs/ds.*patchSize)/2)+0.5)*2)-1; 
        smallBlurPatches = reshape(smallBlurPatches', [patchSize, nPatches]); 
        smallBlurPatches = volresizeSimple(smallBlurPatches, [subsampledSize, nPatches], 'simplelinear');
        smallBlurPatches = reshape(smallBlurPatches, [prod(subsampledSize), nPatches])'; 
        
        % cluster the downsized patches
        smallBlurPatches = bsxfun(@minus, smallBlurPatches, mean(smallBlurPatches, 2));
        gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', gmmTolerance);
        % gmmClust = fitgmdist(smallBlurPatches, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        gmmClust = gmdistribution.fit(smallBlurPatches, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        postVal = gmmClust.posterior(smallBlurPatches);
        [~, clusterIdx] = max(postVal, [], 2);
        
        smallBlurPatches = []; % clear variable
        gmmClust = []; 
        
        save(clusterIdxMat, 'clusterIdx', 'postVal');
    else
        load(clusterIdxMat);
    end

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
    dsPatches = dsPatches(allSamp,:);ls
    wtPatches = wtPatches(allSamp,:); 
    postVal = postVal(allSamp,:); 
    clusterIdx = clusterIdx(allSamp); 
    
    % ppca
    nll = zeros(1, gmmK);
    for k = 1:gmmK
        X0 = dsPatches(clusterIdx==k,:);
        assert(~subtractMean, 'subtractMean not well tested yet');
        if subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        W0 = wtPatches(clusterIdx==k,:); 
        X0(W0<threshold) = nan; 
        
        
        X0 = dsPatches(clusterIdx==k,:); 
        pk = params.ppcaK;
        d = size(X0, 2);
        tic;
        
%         compareCosts(X0, W0, threshold, pk);
        
        % compare init costs.
        if params.ecminit
            warning('detected ecminit');
            X0(W0<threshold) = nan; 
            [~, c]  = ecmninitx(X0, 'twostage');
            [u, s, v] = svd(c); 
            
        elseif isfield(params, 'hackinit') && params.hackinit
            warning('detected hackinit');
            % hacky init: weighted sigma
            m = wmean(X0, W0);
            X0mm = bsxfun(@minus, X0, m);
            W0 = max(W0, 0.00001);
            sq = 0;
            dx = eps;
            for i = 1:size(X0, 1)
                xw = W0(i, :)' .* X0mm(i, :)';
                sq = sq + xw * xw';
                dx = dx + W0(i, :)' * W0(i, :);
            end
            c = sq ./ dx;
            [u, s, v] = svd(c); 
            X0(W0<threshold) = nan; 
            
        else
            [u, s, v] = svd(cov(X0)); 
            X0(W0<threshold) = nan; 
        end
        
        vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
        Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
        opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'final', 'MaxIter', params.maxPPCAiter);
        [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit);
        if subtractMean
            srecon = bsxfun(@minus, S.Recon, mean(S.Recon, 2));
        else
            srecon = S.Recon;
        end
        means(k,:) = mean(srecon);
        sigmas(:,:,k) = cov(srecon);
        fprintf('ppca cluster %d done in %5.3fs (%d patches)\n', k, toc, size(X0, 1));
       
        nll(k) = S.nloglk_new;
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
    nll

    %% save wgmm
    wg = wgmm(means, sigmas, pis);
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params' ); 

    warning('remember to add small diagonal component to sigmas in next code'); 
end

function compareCosts(Y, W, threshold, pk)
    d = size(Y, 2);
    badidx = sum(W, 2) == 0;
    Y(badidx, :) = [];
    W(badidx, :) = [];
    
    [u, s, v] = svd(cov(Y)); 
    vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
    Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
    
    Y(W<threshold) = nan; 
    nlogl1 = ppcax_incomplete_nlogl(Y', nanmean(Y), Winit, vinit);
    
    
    
    
    [mu, c]  = ecmninitx(Y, 'twostage');
    [u, s, v] = svd(c); 
    vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
    Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
    nlogl2 = ppcax_incomplete_nlogl(Y', mu, Winit, vinit);
    
    fprintf('Cov: %f. ECMinit: %f\n', nlogl1, nlogl2);
end




    