function subvol2diagWLSwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
% this file structure is a relic from MICCAI 2016 / Thesis 2016 attempt. We'll keep it for now...
% we're hard-coding dsidx2S init.
% 
% dsSubvolMat - matfile name of subvolume
% wtSubvolMat - matfile name of weights
% iniFilename - ini filename with parameters (input)
% clusterIdxMat - matfile name of cluster assignments (if don't have, where to save)
% wgmmMat - matfile name of output wgmm (output)

    params = ini2struct(iniFilename); 
    
    % rename the parameters
    K = params.nClust;
    patchSize = params.patchSize; 
    nPatches = params.nPatches; 
    
    % volume processing
    volsKeep = params.volsKeep; 
    diffPad = params.volPad;  
    
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
    assert(all(size(dsSubvols) == size(wtSubvols)));
    fprintf('took %5.3f to load the subvolumes\n', toc);

    tic
    % prune volumes
    if isnumeric(volsKeep)
        nSubj = size(dsSubvols, 4);
        nSubjkeep = round(volsKeep .* nSubj);
        fprintf('keeping %d of %d subjects\n', nSubjkeep, nSubj);
        vidx = randsample(nSubj, nSubjkeep);
        
    else % assume it's a file with numbers
        fid = fopen(volsKeep);
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
    
    % crop ds volumes -- we don't want to learn from the edges of the volumes;
    dsSubvolsCrop = cropVolume(dsSubvols, [diffPad+1 1], [volSize-diffPad nSubj]);
    fprintf('took %5.3f to process the %d subvolumes\n', toc, nSubj);
       

    
    %% get patches and init.
    dsPatches = robustVols2lib(dsSubvolsCrop, patchSize);
    whos dsSubvolsCrop
    whos dsPatches
    clear dsSubvolsCrop;

    % Need to select which patches to work with.
    % here we'll use clusterIdxMat as the file to dump the trainsetIdx matrix
    if ~sys.isfile(clusterIdxMat)
        trainsetIdx = randsample(size(dsPatches, 1), nPatches);
        Y = dsPatches(trainsetIdx, :);
        dso = params.dsGmm;

        % this is unnecessary for diag method. we should get rid of it...
        gmmopt = statset('Display', 'iter', 'MaxIter', dso.maxIter, 'TolFun', dso.tol);
        % gmmClust = fitgmdist(smallBlurPatches, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        % gmdist = gmdistribution.fit(Y, gmmK, regstring, 1e-4, 'replicates', 3, 'Options', gmmopt);
        % [~, wgDs] = fitgmdist2wgmmLS(Y, gmmK, 1e-4, 3, gmmopt, dLow);
        [wgDs, wgDsLs] = fitgmdist2wgmmLS(Y, params.nClust, dso.regVal, dso.reps, gmmopt, params.wgmm.dLow);
        
        save(clusterIdxMat, 'trainsetIdx', 'wgDs', 'wgDsLs');
        fprintf('took %5.3f to prepare wgDs\n', toc);
    else
        q = load(clusterIdxMat, 'trainsetIdx', 'wgDs', 'wgDsLs');
        trainsetIdx = q.trainsetIdx;
        wgDs = q.wgDs;
        wgDsLs = q.wgDsLs;
        assert(numel(trainsetIdx) == nPatches, 'The saved number of patches is incorrect');
        
        Y = dsPatches(trainsetIdx, :);
        fprintf('took %5.3f to load wgDs and Idx from %s\n', toc, clusterIdxMat);
    end
    clear dsPatches
    disp(size(Y))
    
    % process weights (later to avoid keeping dsPatches and wtPatches in memory at same time)
    wtSubvolsCrop = cropVolume(wtSubvols, [diffPad+1 1], [volSize-diffPad nSubj]);
    wtPatches = robustVols2lib(wtSubvolsCrop, patchSize);
    sizeWtSubvolsCrop = size(wtSubvolsCrop);
    clear wtSubvolsCrop;
    w = wtPatches(trainsetIdx, :);
    clear wtPatches
    
    

    %% entropy to update W
    
    % get entropy of ds subvolumes
    tic
    enSubvol = dsSubvols*nan; 
    for i = 1:size(enSubvol, 4)
        enSubvol(:,:,:,i) = entropyfilt(dsSubvols(:,:,:,i), getnhood(strel('sphere', 2))); 
    end
    croppedEnSubvols = cropVolume(enSubvol, [diffPad + 1, 1], [volSize - diffPad, nSubj]);
    assert(all(size(croppedEnSubvols) == sizeWtSubvolsCrop));

    % get patches
    enPatchCol = robustVols2lib(croppedEnSubvols, patchSize);
    enDs = enPatchCol(trainsetIdx, :);
    clear enPatchCol;

    % adni
    polyx = [params.entropy.lowEntropy, params.entropy.highEntropy];
    polyy = [params.entropy.lowThr, params.entropy.highThr];
    enFit = polyfit(polyx, polyy, 1);
    wtThrEnDs = within([0.01, params.entropy.highThr], polyval(enFit, enDs));
    
    wtEnDs = w > wtThrEnDs;
    fprintf('mean entropy-based wt: %3.2f\n', mean(wtEnDs(:)))

    data = struct('Y', Y, 'W', double(wtEnDs), 'K', K);
    fprintf('took %5.3f to prepare entropy weighting\n', toc);
    
    %% run wgmm
    wgmmOpts = params.wgmm;
    
    % prepare the wgDsDiag
    wgDsLs.expect = wgDs.expect;
    
    dHigh = size(Y, 2);
    wgDsLsDiag = wgmm(wgDsLs.opts, wgDsLs.params);
    wgDsLsDiag.params.sigma = repmat(eye(dHigh), [1, 1, K]); 
    wgDsLsDiag.params.W = wgDsLsDiag.params.sigma(:, 1:wgmmOpts.dLow, :);
    wgDsLsDiag.expect = wgDs.expect;
    
    % ecm
    wg = wgmmfit(data, 'modelName', 'latentSubspace', 'modelArgs', struct('dopca', wgmmOpts.dLow), ...
        'minIter', wgmmOpts.minIters, 'maxIter', wgmmOpts.maxIters, 'TolFun', wgmmOpts.tolFun, ...
        'verbose', 2, 'replicates', wgmmOpts.reps, ...
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


function lib = robustVols2lib(vols, patchSize)
    try
        lib = patchlib.vol2lib(vols, [patchSize 1]);
    catch err
        fprintf(2, 'vol2lib caught %s\n. Going volume by volume', err.message)
        nVols = size(vols, 4);
        libCell = cell(nVols, 1);
        for i = 1:nVols
            libCell{i} = patchlib.vol2lib(vols(:,:,:,i), patchSize);
        end
        lib = cat(1, libCell{:});
    end
end

