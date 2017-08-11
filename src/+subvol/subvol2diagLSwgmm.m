function subvol2diagLSwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
% this file structure is a relic from MICCAI 2016 / Thesis 2016 attempt. We'll keep it for now...
% we're hard-coding dsidx2S init.
% 
% dsSubvolMat - matfile name of subvolume
% wtSubvolMat - matfile name of weights
% iniFilename - ini filename with parameters (input)
% clusterIdxMat - matfile name of cluster assignments (if don't have, where to save)
% wgmmMat - matfile name of output wgmm (output)

    params = ini2struct(iniFilename); 
    
    percVolsKeep = params.percVolsKeep; 

    % rename the parameters
    gmmK = params.gmmK; 
    patchSize = params.patchSize; 
    nPatches = params.nPatches; 
    
    gmmTolerance = params.gmmTolerance; 
    regstring = params.regstring;
    
    diffPad = params.volPad;
    dLow = params.dLow;
    maxIters = params.maxIters;
    minIters = params.minIters;
    TolFun = params.TolFun;
    wthr = params.threshold;
    
    growMethod = params.growMethod;
    
    nDims = numel(patchSize);

    %% load subvolumes
    tic
    if ischar(dsSubvolMat)
        q = load(dsSubvolMat, 'subvolume'); 
        dsSubvols = q.subvolume;
    else
        dsSubvols = dsSubvolMat;
    end
    if ischar(wtSubvolMat)
        q = load(wtSubvolMat, 'subvolume'); 
        wtSubvols = q.subvolume;
    else 
        wtSubvols = wtSubvolMat;
    end
    clear q;
    fprintf('took %5.3f to load the subvolumes\n', toc);

    tic
    % prune volumes
    if isnumeric(percVolsKeep)
        nSubj = size(dsSubvols, 4);
        nSubjkeep = round(percVolsKeep .* nSubj);
        vidx = randsample(nSubj, nSubjkeep);
    else % assume it's a file with numbers
        fid = fopen(percVolsKeep);
        C = textscan(fid, '%d');
        fclose(fid);
        vidx = cat(1, C{:});
        nSubjkeep = numel(vidx);
    end
    dsSubvols = dsSubvols(:,:,:,vidx);
    wtSubvols = wtSubvols(:,:,:,vidx);
    nSubj = nSubjkeep;
    volSize = size(dsSubvols);
    volSize = volSize(1:3);
    
    if isscalar(diffPad)
        diffPad = diffPad * ones(1, 3);
    end
    
    % crop volumes -- we don't want to learn from the edges of the volumes;
    dsSubvolsCrop = cropVolume(dsSubvols, [diffPad+1 1], [volSize-diffPad nSubj]);
    wtSubvolsCrop = cropVolume(wtSubvols, [diffPad+1 1], [volSize-diffPad nSubj]);
    
    fprintf('took %5.3f to process the %d subvolumes\n', toc, nSubj);
    
    %% get patches and init.
    dsPatches = patchlib.vol2lib(dsSubvolsCrop, [patchSize 1]);
    wtPatches = patchlib.vol2lib(wtSubvolsCrop, [patchSize 1]);
    

    
    % Need to select which patches to work with.
    % here we'll use clusterIdxMat as the file to dump the trainsetIdx matrix
    if ~sys.isfile(clusterIdxMat)
        trainsetIdx = randsample(size(dsPatches, 1), nPatches);
        Y = dsPatches(trainsetIdx, :);
        data = struct('Y', Y, 'W', wtPatches(trainsetIdx, :) > wthr, 'K', gmmK);

        
        gmmopt = statset('Display', 'iter', 'MaxIter', 3, 'TolFun', gmmTolerance);
        % gmmClust = fitgmdist(smallBlurPatches, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        % gmdist = gmdistribution.fit(Y, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        % [~, wgDs] = fitgmdist2wgmmLS(Y, gmmK, 1e-4, 3, gmmopt, dLow);
        [wgDs, wgDsLs] = fitgmdist2wgmmLS(Y, gmmK, 1e-4, 3, gmmopt, dLow);
        
        save(clusterIdxMat, 'trainsetIdx', 'wgDs', 'wgDsLs');
    else
        q = load(clusterIdxMat, 'trainsetIdx', 'wgDs', 'wgDsLs');
        trainsetIdx = q.trainsetIdx;
        wgDs = q.wgDs;
        wgDsLs = q.wgDsLs;
        assert(numel(trainsetIdx) == nPatches, 'The saved number of patches is incorrect');
        
        Y = dsPatches(trainsetIdx, :);
        data = struct('Y', Y, 'W', wtPatches(trainsetIdx, :) > wthr, 'K', gmmK);
    end
        
    %% run    
    
        % prepare the wgDsDiag
        wgDsLs.expect = wgDs.expect;
        
        dHigh = size(Y, 2);
        wgDsLsDiag = wgmm(wgDsLs.opts, wgDsLs.params);
        wgDsLsDiag.params.sigma = repmat(eye(dHigh), [1,1,gmmK]); 
        wgDsLsDiag.params.W = wgDsLsDiag.params.sigma(:, 1:dLow, :);
        wgDsLsDiag.expect = wgDs.expect;
    
    % ecm
    itersteps = params.ppcaMinK:params.ppcaKskip:params.ppcaMaxK;
    wg = wgmmfit(data, 'modelName', 'latentSubspace', 'modelArgs', struct('dopca', dLow), ...
        'minIter', minIters, 'maxIter', maxIters, 'TolFun', TolFun, 'verbose', 2, 'replicates', 1, ...
        'init', 'wgmm', 'initArgs', struct('wgmm', wgDsLsDiag));

    wgwhole = wg;
    wg = wgmm(wg.opts, wg.params);
    if isfield(wg.params, 'sigma')
        wg.params = rmfield(wg.params, 'sigma');
    end
    
    %% save
    if params.saveLargeWg
        save(wgmmMat, 'wg', 'wgwhole', 'trainsetIdx', '-v7.3'); 
    else
        save(wgmmMat, 'wg', 'trainsetIdx', '-v7.3'); 
    end
end
