function wg = subvol2iterppcatest(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename, optionalSubvolSaveMatFileName)
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
    
    

    wparams = struct();
    for k = 1:gmmK
            assert(~subtractMean, 'subtractMean not well tested yet');

            X0 = dsPatches(clusterIdx==k,:); 
            W0 = wtPatches(clusterIdx==k,:); 
            %params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK
            d = size(X0, 2);
            c = cov(X0); % initial covariance

            wparams.mu(k,:) = mean(X0);
            wparams.sigma(:,:,k) = c;
            
%             pk = 1;
%             [u, s, ~] = svd(c); 
%             vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
%             Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
%             
%             wparams.sigmasq(k) = vinit;
%             wparams.W(:,:,k) = Winit;
    end
    wparams.pi = sum(postVal) ./ sum(postVal(:));
    
    % test
    mi = argmax(postVal, [], 2);
    piv = postVal * 0;
    idx = sub2ind(size(postVal), (1:size(postVal, 1))', mi(:));
    piv(idx) = 1;
    data = struct('Y', dsPatches, 'W', wtPatches>=threshold, 'K', gmmK);
    wginit = wgmm(wgmm.optionDefaults, wparams);
    wginit.expect.gammank = piv;
    qwg = wgmmfit(data, 'replicates', 1, 'modelName', 'latentSubspace', ...
        'modelArgs', struct('dopca', inf), 'maxIter', params.maxPPCAiter, 'minIter', params.maxPPCAiter, 'TolFun', 0.00001, ...
        'verbose', 2, 'init', 'latentSubspace-iterds-randv', ...
        'initargs', struct('wgmm', wginit, 'dopcas', params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK,...
        'growW', 'recompute'));
    save('qwg', 'qwg');
    wg = qwg;
    
    % ppca
    dsadsadas
    for k = 1:gmmK
        assert(~subtractMean, 'subtractMean not well tested yet');
        
        X0 = dsPatches(clusterIdx==k,:); 
        W0 = wtPatches(clusterIdx==k,:); 
        %params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK
        d = size(X0, 2);
        c = cov(X0); % initial covariance
        
        means0(k,:) = mean(X0);
        means(k,:) = means0(k,:);
        sigmas0(:,:,k) = c;      
        
        X0(W0<threshold) = nan; 
        tic;
        for pk = params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK
            % compute initialization
            [u, s, v] = svd(c); 
            vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
            Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
            opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'final', 'MaxIter', params.maxPPCAiter);
            %[COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit); %, 'v0', vinit); % v0 added later!!!
            
            
            % for comparison, let's set a common v0
            v0init = rand;
            [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit, 'v0', v0init, 'mu0', means(k,:)'); % v0 added later!!!
            
            % test wgmm
            data = struct('Y', X0, 'W', W0>=threshold, 'K', 1);
            wginit = wgmm(wgmm.optionDefaults, struct('mu', means(k,:), 'W', Winit, 'sigmasq', v0init, 'pi', 1));
            qwglocal = wgmmfit(data, 'replicates', 1, 'modelName', 'latentSubspace', ...
                'modelArgs', struct('dopca', pk), 'maxIter', params.maxPPCAiter, 'TolFun', 0.00001, ...
                'verbose', 2, 'init', 'wgmm', 'initargs', struct('wgmm', wginit));
% %             
            % reconstruction
            srecon = S.Recon;
            means(k,:) = mean(srecon);
            means(k,:) = MU(:)';
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

    %% save wgmm
    wg = wgmm([], struct('mu', means, 'sigma', sigmas, 'pi', pis));
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params', 'Ssave'); 

    warning('remember to add small diagonal component to sigmas in next code'); 
