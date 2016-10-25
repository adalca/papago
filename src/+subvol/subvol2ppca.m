function subvol2ppca(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
% dsSubvolMat - matfile name of subvolume
% wtSubvolMat - matfile name of weights
% iniFilename - ini filename with parameters (input)
% clusterIdxMat - matfile name of cluster assignments (if don't have, where to save)
% wgmmMat - matfile name of output wgmm (output)

    params = ini2struct(iniFilename); 

    % rename the parameters
    nDims = numel(params.maxPatchSize);

    % load subvolumes
    tic
    q = load(dsSubvolMat, 'subvolume'); 
    dsSubvols = q.subvolume;
    q = load(wtSubvolMat, 'subvolume'); 
    wtSubvols = q.subvolume;
    clear q;
    fprintf('took %5.3f to load the subvolumes\n', toc);
    
    % compute some important volume statistics
    volSize = size(dsSubvols);
    nSubj = volSize(nDims + 1);
    
    for sc = 1:params.ds
        patchSize = round((params.maxPatchSize - params.minPatchSize) .* sc ./ params.ds + params.minPatchSize);
        params.patchSize = patchSize;
        params.superPatchSize  = params.patchSize + 4;
        
        diffPad = (params.superPatchSize-patchSize)/2;
        assert(any(isIntegerValue(diffPad)), 'difference between superPatchSize and patchSize should be even');
        params.diffPad = diffPad; params.nSubj = nSubj;

        % prepare gmm index (blur subvolume and train gmm)
        % TODO: only do this @ first level? from then on do gmm at previous level (?).
        [clusterIdx, postVal, pis] = initgmm(dsSubvols, wtSubvols, clusterIdxMat, params);

        % crop subvolumes and get patches of size patchSize (e.g. 9x9x9)
        dsPatches = cropVolume(dsSubvols, [diffPad+1 1], [size(dsSubvols(:,:,:,1))-diffPad nSubj]);
        [dsPatches, ~, ~, srcgridsize] = patchlib.vol2lib(dsPatches, [patchSize 1]);
        wtPatches = cropVolume(wtSubvols, [diffPad+1 1], [size(wtSubvols(:,:,:,1))-diffPad nSubj]);
        wtPatches = patchlib.vol2lib(wtPatches, [patchSize 1]);
        assert(numel(clusterIdx) == size(dsPatches, 1));

        % sample and initlize covariance/mean
        [means0, sigmas0, allSamp] = init(dsPatches, wtPatches, params);
        dsPatches = dsPatches(allSamp,:); 
        wtPatches = wtPatches(allSamp,:); 
        postVal = postVal(allSamp,:); 
        clusterIdx = clusterIdx(allSamp); 

        % compute parameters
        [means, sigmas] = mstep(dsPatches, wtPatches, means0, sigmas0, clusterIdx, params);

        % estimate values for all missing patches.
        filledPatches = fill(dsPatches, wtPatches, means, sigmas, clusterIdx, params);

        % quilt patches
        quiltedVol = patchlib.quilt(filledPatches, srcgridsize, [patchSize, 1]);
        
        % prepare next scale
        if sc < params.ds
            % get size of next scale
            nextScaleSize = round(size(dsSubvols) * sc ./ params.ds);
            % upsample volume
            scQuiltedVol = volresize(quiltedVol, nextScaleSize);
            % downsample volume from orig to next scale
            vol = subvolsresize(dsSubvols, wtSubvols, nextScaleSize);
            
            % combine volumes
            scQuiltedVol
        end
    end
    
    
    % save final wgmm
    wg = wgmm(means, sigmas, pis);
    save(wgmmMat, 'wg', 'allSamp', 'clusterIdx', 'postVal', 'params' ); 

    warning('remember to add small diagonal component to sigmas in next code'); 

end


function [filledPatches] = fill(dsPatches, wtPatches, means, sigmas, clusterIdx, params)
% perform mstep (ecm or ppca) at a g

    gmmK = params.gmmK; 
    filledPatches = zeros(size(dsPatches));
    
    for k = 1:gmmK

        % get blurry patches
        X0 = dsPatches(clusterIdx == k, :);
        if params.subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        % get weights
        W0 = wtPatches(clusterIdx == k, :);
        % set up downsampled data with nans in any region with W < threshold
        X0(W0<threshold) = nan;

        % fill in patches with p(a|b)
        filledPatches(clusterIdx == k, :) = inpaintWithGaussConditional(X0, means(k,:), sigmas(:,:,k));
    end
end

function [means, sigmas] = mstep(dsPatches, wtPatches, means0, sigmas0, clusterIdx, params)
% perform mstep (ecm or ppca) at a g

    gmmK = params.gmmK; 
    patchSize = params.patchSize; 
    
    means = zeros(gmmK, prod(patchSize)); 
    sigmas = zeros(prod(patchSize), prod(patchSize), gmmK); 
    for k = 1:gmmK

        % get blurry patches
        X0 = dsPatches(clusterIdx == k, :);
        if params.subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        % get weights
        W0 = wtPatches(clusterIdx == k, :);
        % set up downsampled data with nans in any region with W < threshold
        X0(W0<threshold) = nan;

        % if less than hackNum non NaN elements in column then add in some new points
        [~, idxWvals] = sort(W0, 'descend');
        inds = sub2ind( size(W0), idxWvals(1:params.minNonNan,:), ndgrid(1:size(W0,2),1:params.minNonNan)' ); 
        X0( inds(:) ) = X0( inds(:) );

        % run ecm 
        if strcmp(params.method, 'ecm')
            tic;
            [means(k,:), sigmas(:,:,k), Zout] = ecmnmlex(X0, 'twostage', params.maxECMIter, params.tolerance, means0(k,:), sigmas0(:,:,k));
            fprintf('took %5.3f for ecm cluster %d\n', toc, k);

%             if ~params.subtractMean
%                 warning('hacky computation of non-subtracted-stats for now: Basically, since all other components assume extracted mean, we''re computing the stats based on the infered data. This is very very ugly :(');
%                 means(k,:) = mean(bsxfun(@minus, Zout, mean(Zout, 2))); 
%                 sigmas(:,:,k) = cov(bsxfun(@minus, Zout, mean(Zout, 2))) + eye(size(X0,2)) * 0.00001;
%             end
%             assert(isclean(sigmas));
            
        else, assert(strcmp(params.method, 'ppca'))
             assert(~params.subtractMean, 'ppca has not been tested with subtractMean on');
    
        
            pk = params.ppcaK;
            d = size(X0, 2);
            tic;
            [u, s, v] = svd(cov(X0)); 

            X0(W0<threshold) = nan; 
            vinit = (1 ./ (d-pk)) * sum(diag(s((pk+1):d, (pk+1):d)));
            Winit = u(:, 1:pk) * (sqrt(s(1:pk, 1:pk)) - vinit * eye(pk));
            opts = struct('TolFun', params.ppcaTolFun, 'TolX', params.ppcaTolX, 'Display', 'iter', 'MaxIter', params.maxPPCAiter);
            [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, pk, 'Options', opts, 'W0', Winit);
%             if params.subtractMean
%                 srecon = bsxfun(@minus, S.Recon, mean(S.Recon, 2));
%             else
                srecon = S.Recon;
%             end
            means(k,:) = mean(srecon);
            sigmas(:,:,k) = cov(srecon) + eye(size(X0,2)) * 0.00001;
            fprintf('ppca cluster %d done in %5.3fs\n', k, toc);
        end
    end
end


function [means0, sigmas0, allSamp] = init(dsPatches, wtPatches, params)

    gmmK = params.gmmK; 
    patchSize = params.patchSize; 

    means0 = zeros(gmmK, prod(patchSize)); 
    sigmas0 = zeros(prod(patchSize), prod(patchSize), gmmK);
    
    nPatches = size(dsPatches,1); 
    randomIdx = randperm(nPatches);
    clusterIdxPerm = clusterIdx(randomIdx); 
    allSamp = []; 
    
    for k = 1:gmmK
        X0 = dsPatches(clusterIdx==k,:); 
        if params.subtractMean
            X0 = bsxfun(@minus, X0, mean(X0, 2));
        end
        W0 = wtPatches(clusterIdx==k,:); 
        X0(W0<threshold) = nan; 
        
        tic;
        warning('using non mean-subtracted');
        [means0(k,:), sigmas0(:,:,k)] = ecmninitx(X0, 'twostage'); 
        fprintf('took %5.3f to init ecm\n', toc);
        
        samp = randomIdx(clusterIdxPerm==k); 
        allSamp = [allSamp samp(1:min(params.nPatchesThresh, length(samp)))]; 
    end
end

function [clusterIdx, postVal, pis] = initgmm(dsSubvols, wtSubvols, clusterIdxMat, params)

    % gmm params
    gmmK = params.gmmK; 
    ds = params.ds;
    smallUs = params.smallUs; 
    gmmTolerance = params.gmmTolerance;
    RECOMPUTECLUSTERS = params.RECOMPUTECLUSTERS;
    regstring = params.regstring;
    nSubj = params.nSubj;
    diffPad = params.diffPad;

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
end

function [vol, wvol] = subvolsresize(dsSubvols, wtSubvols, newSize)
    sigmaBlur = 1/3 * mean(newSize ./ size(dsSubvols));

    % weight-based blurring
    dsBlurvol = volblur(dsSubvols .* wtSubvols, sigmaBlur); 
    wtBlurvol = volblur(wtSubvols, sigmaBlur); 
    vol = volresize(dsBlurvol./wtBlurvol, newSize);  
    wvol = volresize(wtBlurvol, newSize);  
    warning('perhaps todo: rescale weight?');
end
