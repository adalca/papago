function subvol2ecmwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
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
    minNonNan = params.minNonNan; 
    nPatchesThresh = params.nPatchesThresh; 
    maxECMIter = params.maxECMIter; 
    tolerance = params.tolerance; 
    RECOMPUTECLUSTERS = params.RECOMPUTECLUSTERS; 

    %% load subvolumes
    q = load(dsSubvolMat); 
    dsSubvols = q.subvolumes;
    q = load(wtSubvolMat); 
    wtSubvols = q.subvolumes;

    %% split subvolumes into super patches (e.g. 13x13x13)
    dsPatches = patchlib.vol2lib(dsSubvols, [superPatchSize 1]); 
    wtPatches = patchlib.vol2lib(wtSubvols, [superPatchSize 1]);

    if (~exist(clusterIdxMat, 'file') || RECOMPUTECLUSTERS)
        % train clusters

        % downsize the patches using the weights
        warning('only using the sizes of the first dimension'); 
        [dowsizePatches] = appxDownsizePatches(dsPatches, wtPatches, superPatchSize(1), patchSize(1), ds, smallUs);

        % cluster the downsized patches
        meanAdjDownsizePatches = bsxfun(@minus, dowsizePatches, mean(dowsizePatches, 2));
        gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', tolerance);
        gmmClust = fitgmdist(meanAdjDownsizePatches, gmmK, 'regularizationValue', 1e-4, 'replicates', 3, 'Options', gmmopt);
        postVal = gmmClust.posterior(meanAdjDownsizePatches);
        [~, clusterIdx] = max(postVal, [], 2);

        save(clusterIdxMat, 'clusterIdx', 'postVal');

    else
        load(clusterIdxMat);
    end


    %% crop patches to 9x9x9
    % reshape to be volumes
    N = size(dsPatches,1); 
    dsPatchesVols = reshape(dsPatches, [superPatchSize N]);
    wtPatchesVols = reshape(wtPatches, [superPatchSize N]);

    % crop out the cetner of the volume
    diffPad = (superPatchSize-patchSize)/2;
    dsPatchesVolsCrop = cropVolume(dsPatchesVols, [diffPad+1 1], [superPatchSize-diffPad N]); 
    wtPatchesVolsCrop = cropVolume(wtPatchesVols, [diffPad+1 1], [superPatchSize-diffPad N]); 

    % reshape back
    dsPatchesVolsCrop = reshape(dsPatchesVolsCrop, [prod(patchSize) N])'; 
    wtPatchesVolsCrop = reshape(wtPatchesVolsCrop, [prod(patchSize) N])'; 


    %% run ecm
    means = zeros(gmmK, prod(patchSize)); 
    sigmas = zeros(prod(patchSize), prod(patchSize), gmmK); 

    for k = 1:gmmK

        % get blurry patches
        X0 = dsPatchesVolsCrop(clusterIdx == k, :);
        X0 = bsxfun(@minus, X0, mean(X0, 2));

        % get weights
        W = wtPatchesVolsCrop(clusterIdx == k, :);

        % set up downsampled data with nans in any region with W < threshold
        X0nans = X0;
        X0nans(W<threshold) = nan;

        % get initial mean and covariance
        [Mean0, Covar0] = ecmninit(X0nans, 'twostage'); 

        % sample patches 
        nclustpatches = sum(clusterIdx==k);
        samp = randperm(nclustpatches,min(nPatchesThresh,nclustpatches));
        X0nansampled = X0nans(samp,:);
        X0sampled = X0(samp,:); 

        % if less than hackNum non NaN elements in column then add in some new points
        [~, idxWvals] = sort(W, 'descend');
        nanColumns = sum(~isnan(X0nansampled),1) < minNonNan;
        inds = sub2ind( size(X0nansampled),  idxWvals(1:minNonNan, nanColumns), repmat( find(nanColumns == 1), [minNonNan, 1]) );
        X0nansampled( inds(:) ) = X0sampled( inds(:) );

        % run ecm 
        [means(k,:), sigmas(:,:,k)] = ecmnmlex(X0nansampled, 'twostage', maxECMIter, tolerance, Mean0, Covar0);

    end

    %% output wgmm
    wg = wgmm(means, sigmas, hist(clusterIdx, 1:gmmK)./numel(clusterIdx));
    save(wgmmMat, 'wg', 'clusterIdx', 'postVal', 'params' ); 

    warning('remember to add small diagonal component to sigmas in next code'); 

