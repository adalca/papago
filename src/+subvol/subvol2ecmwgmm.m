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
    
    means0 = zeros(gmmK, prod(patchSize)); 
    sigmas0 = zeros(prod(patchSize), prod(patchSize), gmmK);
    
    for k = 1:gmmK
        X0 = dsPatches(clusterIdx==k,:); 
        if subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        W0 = wtPatches(clusterIdx==k,:); 
        X0(W0<threshold) = nan; 
        
        if params.ppcaInit
            X0 = dsPatches(clusterIdx==k,:); 
            
            
            pk = params.ppcaK;
            d = size(X0, 2);
            tic;
            [u, s, v] = svd(cov(X0)); 
            
            X0(W0<threshold) = nan; 
            vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
            Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
            opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'iter', 'MaxIter', params.maxPPCAiter);
            [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit);
            if subtractMean
                srecon = bsxfun(@minus, S.Recon, mean(S.Recon, 2));
            else
                srecon = S.Recon;
            end
            means0(k,:) = mean(srecon);
            sigmas0(:,:,k) = cov(srecon) + eye(size(X0,2)) * 0.00001;
            fprintf('ppca cluster %d done in %5.3fs\n', k, toc);
        else
            tic;
            warning('using non mean-subtracted');
            [means0(k,:), sigmas0(:,:,k)] = ecmninitx(X0, 'twostage'); 
            fprintf('took %5.3f to init ecm\n', toc);
        end
    end
    
    %% run ecm 
    means = zeros(gmmK, prod(patchSize)); 
    sigmas = zeros(prod(patchSize), prod(patchSize), gmmK); 
    for k = 1:gmmK

        % get blurry patches
        X0 = dsPatches(clusterIdx == k, :);
        if subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end

        % get weights
        W0 = wtPatches(clusterIdx == k, :);

        % set up downsampled data with nans in any region with W < threshold
        X0 = X0;
        X0(W0<threshold) = nan;

        % if less than hackNum non NaN elements in column then add in some new points
        [~, idxWvals] = sort(W0, 'descend');
        inds = sub2ind( size(W0), idxWvals(1:minNonNan,:), ndgrid(1:size(W0,2),1:minNonNan)' ); 
        X0( inds(:) ) = X0( inds(:) );

        % run ecm 
        tic;
        [means(k,:), sigmas(:,:,k), Zout] = ecmnmlex(X0, 'twostage', maxECMIter, tolerance, means0(k,:), sigmas0(:,:,k));
        fprintf('took %5.3f for ecm cluster %d\n', toc, k);
        
%         if ~subtractMean
%             warning('hacky computation of non-subtracted-stats for now: Basically, since all other components assume extracted mean, we''re computing the stats based on the infered data. This is very very ugly :(');
%             means(k,:) = mean(bsxfun(@minus, Zout, mean(Zout, 2))); 
%             sigmas(:,:,k) = cov(bsxfun(@minus, Zout, mean(Zout, 2)));
%         end
        assert(isclean(sigmas));
    end

    %% save wgmm
    wg = wgmm(means, sigmas, pis);
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params' ); 

    warning('remember to add small diagonal component to sigmas in next code'); 


