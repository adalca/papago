function subvol2iterppcawgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename, optionalSubvolSaveMatFileName)
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
    dsPatcheso = cropVolume(dsSubvols, [diffPad+1 1], [size(dsSubvols(:,:,:,1))-diffPad nSubj]);
    [dsPatcheso, ~, ~, srcgridsize] = patchlib.vol2lib(dsPatcheso, [patchSize 1]);
    
    wtPatcheso = cropVolume(wtSubvols, [diffPad+1 1], [size(wtSubvols(:,:,:,1))-diffPad nSubj]);
    wtPatcheso = patchlib.vol2lib(wtPatcheso, [patchSize 1]);
    
    
    %% sample and initlize covariance/mean

    means = zeros(gmmK, prod(patchSize)); 
    sigmas = zeros(prod(patchSize), prod(patchSize), gmmK);
    
    % sample first.
    nPatches = size(dsPatcheso,1); 
    randomIdx = randperm(nPatches);
    clusterIdxPerm = clusterIdx(randomIdx); 
    allSamp = []; 
    for k = 1:gmmK
        samp = randomIdx(clusterIdxPerm==k); 
        locsamp = samp(1:min(nPatchesThresh,length(samp)));
        allSamp = [allSamp locsamp]; 
    end
    dsPatches = dsPatcheso(allSamp,:); 
    wtPatches = wtPatcheso(allSamp,:); 
    postVal = postVal(allSamp,:); 
    clusterIdx = clusterIdx(allSamp); 
    
    % ppca
    for k = 1:gmmK
        assert(~subtractMean, 'subtractMean not well tested yet');
        
        X0 = dsPatches(clusterIdx==k,:); 
        W0 = wtPatches(clusterIdx==k,:); 
        params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK
        d = size(X0, 2);
        c = cov(X0); % initial covariance
        
        means0(k,:) = mean(X0);
        sigmas0(:,:,k) = c;
        
        
        X0(W0<threshold) = nan; 
        tic;
        for pk = params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK
            % compute initialization
            [u, s, v] = svd(c); 
            vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
            Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
            opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'final', 'MaxIter', params.maxPPCAiter);
            [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit);
            1
            % reconstruction
            srecon = S.Recon;
            means(k,:) = mean(srecon);
            sigmas(:,:,k) = cov(srecon);
            sigmas(:,:,k) = S.W * S.W' + S.v * eye(size(S.W,1));
            c = sigmas(:,:,k);
        end
        Ssave{k} = S;
        
%         if ~subtractMean
%             % if subtractMean is off, then we did not subtract the 
%             % patchwise mean to create these covariances. However, the
%             % reconstruction code right now assumes we do, and expects mean
%             % and cov as such
%             warning('hacky subtractMean fix. This is very very ugly :(');
%             means(k,:) = mean(bsxfun(@minus, srecon, mean(srecon, 2))); 
%             sigmas(:,:,k) = cov(bsxfun(@minus, srecon, mean(srecon, 2)));
%         end
%         
        fprintf('iter-ppca cluster %d done in %5.3fs\n', k, toc);
    end
    
    %% recon in atl space
    if exist('optionalSubvolSaveMatFileName', 'var')
        tic;
        Y_hat = zeros(size(dsPatcheso));
        for k = 1:gmmK 
            X0 = dsPatcheso(clusterIdx==k,:); 
            W0 = wtPatcheso(clusterIdx==k,:); 

            X0(W0<threshold) = nan; 
            mu = mean(Ssave{k}.Recon);
            Y_hat(clusterIdx == k, :) = ppcax_recon(X0, Ssave{k}.W, Ssave{k}.v, mu);
        end

        for i = 1:srcgridsize(end)
            nVols = size(dsSubvols, 4);
            nVoxels = size(dsPatcheso, 1) ./ nVols;
            idx = ((i-1)* nVoxels + 1):i*nVoxels;
            reconVols(:,:,:,i) = patchlib.quilt(Y_hat(idx, :), [srcgridsize(1:3), 1], [patchSize, 1]);
            i
        end
        
        reconWeight = patchlib.quilt(double(Y_hat(idx, :) *0 + 1), [srcgridsize(1:3), 1], [patchSize, 1]);
        
        reconMasks = wtPatcheso;
        q = load(dsSubvolMat, 'nfo');
        reconLoc = q.nfo.subvolLoc;
        save(optionalSubvolSaveMatFileName, 'reconVols', 'reconMasks', 'reconLoc', 'reconWeight', '-v7.3');
        toc
    end
    
    %% subject-specific manual testing
    doingmanualtesting = false;
    if doingmanualtesting 
        Y_hat = zeros(size(dsPatcheso));
        for k = 1:gmmK 
            X0 = dsPatcheso(clusterIdx==k,:); 
            W0 = wtPatcheso(clusterIdx==k,:); 

            X0(W0<threshold) = nan; 
            mu = mean(Ssave{k}.Recon);
            Y_hat(clusterIdx == k, :) = ppcax_recon(X0, Ssave{k}.W, Ssave{k}.v, mu);
        end

        volno = 23;
        nVols = size(dsSubvols, 4);
        nVoxels = size(dsPatcheso, 1) ./ nVols;

        idx = ((volno-1)* nVoxels + 1):volno*nVoxels;


        volsppca = patchlib.quilt(Y_hat(idx, :), [srcgridsize(1:3), 1], [patchSize, 1]);
        %view3D(vol);

        Y_hat2 = zeros(size(dsPatcheso));
        Y_hat3 = zeros(size(dsPatcheso));
        for i = idx
            k = clusterIdx(i);

            M0 = wtPatcheso(i,:)>=threshold;
            x = dsPatcheso(i,:);
            x(~M0) = nan; 

            Y_hat2(i, :) = paffine.reconSubjPatch(x, logical(M0), logical(M0*0+1), means(k, :)', sigmas(:,:,k)); 
            Y_hat3(i, :) = paffine.reconSubjPatch(x, logical(M0), logical(M0*0+1), means0(k, :)', sigmas0(:,:,k)); 
            i
        end
        volppcarecon = patchlib.quilt(Y_hat2(idx, :), [srcgridsize(1:3), 1], [patchSize, 1]);
        vol0recon = patchlib.quilt(Y_hat3(idx, :), [srcgridsize(1:3), 1], [patchSize, 1]);
        %view3D(vol);

        volds = patchlib.quilt(dsPatcheso(idx,:), [srcgridsize(1:3), 1], [patchSize, 1]);
        volwt = patchlib.quilt(wtPatcheso(idx,:), [srcgridsize(1:3), 1], [patchSize, 1]);

        view3Dopt(volds, volwt, volsppca, volppcarecon, vol0recon);
    end

    %% save wgmm
    wg = wgmm(means, sigmas, pis);
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params', 'Ssave'); 

    warning('remember to add small diagonal component to sigmas in next code'); 
