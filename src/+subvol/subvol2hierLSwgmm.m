function [wgDs, vols] = subvol2hierLSwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
% this file structure is a relic from MICCAI 2016 / Thesis 2016 attempt. We'll keep it for now...
% run hierarchy
%
% see also: testHierarchy
%
% TODO: should really do more of the blur & rotate plan

% should have dsSubvols, wtSubvols, isoSubvols already loaded, as well as interpData

% NOTE: vols returned is from the previous (not final) iteration!.


    %% Parse parameters    
    params = ini2struct(iniFilename); 

    % patch parameters
    scSizes = params.scSizes;               % 5:2:9;
    nPatches = params.nPatches;             % 25000
    
    % volume parameters
    blurWindow = params.blurWindow;         % [15, 15, 15, 1];       % blur window (at main scale)
    percVolsKeep = params.percVolsKeep; 
    diffPad = params.volPad;                % usually 2 (?)
    
    % gmm parameters
    gmmK = params.gmmK; 
    gmmTolerance = params.gmmTolerance;     % 0.001
    regVal = params.regVal;                 % 1e-4;
    wgmmMethod = params.wgmmMethod;         % 'LS_diagw', 'LS_diag', 'LS_ds'
    dLow = params.dLow;                     % usually 21
    maxIters = params.maxIters;             % usually 15-21
    
    % other
    origwthr = params.threshold;            % normally: something like 0.5. LS_diagw: 0.5 for 50, 0 for Full
    reconModel = params.reconModel;         % 'latentMissing', 'latentSubspace'
    
    % compute some useful parameters
    patchSize = [1,1,1] * scSizes(end);
    minPathSize = [1, 1, 1] * scSizes(1);


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

    % prune volumes
    nSubj = size(dsSubvols, 4);
    nSubjkeep = round(percVolsKeep .* nSubj);
    vidx = randsample(nSubj, nSubjkeep);
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

    
    %% parameters

    % percent known data (percent data to trust
    percKnown = mean(wtSubvols(:) > origwthr); % note: here we look at the original weight matrix (original scale.
    % percKnown = 0.05;
    
    % prepare wgmm
    dsgmmopt = statset('Display', 'final', 'MaxIter', 15, 'TolFun', gmmTolerance);
    lsopts = {'modelName', 'latentSubspace', 'modelArgs', struct('dopca', dLow), ...
        'maxIter', maxIters, 'TolFun', gmmTolerance, 'verbose', 2};

    %% Run hierarchy
    
    nSubj = size(dsSubvols, 4);

    % prepare volumes
    scAtlPatchSize = minPathSize;
    blurSigma = 1/3 * patchSize ./ minPathSize;
    tic;
    vols = hierarchyPrepVolumesSimple(dsSubvols, wtSubvols, ...
        blurSigma, blurWindow, diffPad, patchSize, scAtlPatchSize);
    fprintf('Vol extraction took %3.2f\n', toc);

    % loop over scales
    for si = 1:numel(scSizes)
        scAtlPatchSize = [1,1,1]*scSizes(si);

        % extract patches from volumes
        mxPatchesSmall = patchlib.vol2lib(vols.mxSubvolSmallCrop, [scAtlPatchSize, 1]);
        wtPatchesSmall = patchlib.vol2lib(vols.wtSubvolSmallCrop, [scAtlPatchSize, 1]);

        % sample patches for wgmm.
        dataridx = randsample(size(mxPatchesSmall, 1), nPatches);
        wthr = prctile(vols.wtSubvolSmallCrop(:), 100-100*percKnown);
        yds = mxPatchesSmall(dataridx, :);

        switch wgmmMethod
            case 'LS_ds'
                [~, wgDsLs] = fitgmdist2wgmmLS(yds, gmmK, regVal, 3, dsgmmopt, dLow);
                
                Ws = double(wtPatchesSmall(dataridx, :) > wthr);
                wdata = struct('Y', yds, 'W', Ws, 'K', gmmK);
                paramsLS_ds = [lsopts, 'replicates', 1, 'init', 'wgmm', 'initArgs', struct('wgmm', wgDsLs)];
                wgDs = wgmmfit(wdata, paramsLS_ds{:});
            
            case 'LS_diag'
                % run ds wgmm
                [~, wgDsLs] = fitgmdist2wgmmLS(yds, gmmK, regVal, 3, dsgmmopt, dLow);
                % put diagonal along
                wgDsLsDiag = wgmm(wgDsLs.opts, wgDsLs.params);
                wgDsLsDiag.params.sigma = repmat(eye(prod(scAtlPatchSize)), [1,1,gmmK]); 
                wgDsLsDiag.params.W = wgDsLsDiag.params.sigma(:, 1:dLow, :);

                % run LS_diag wgmm on top of ds gmm.
                Ws = double(wtPatchesSmall(dataridx, :) > wthr);
                wdata = struct('Y', yds, 'W', Ws, 'K', gmmK);
                paramsLS_diag = [lsopts, 'init', 'wgmm', 'initArgs', struct('wgmm', wgDsLsDiag)];
                wgDs = wgmmfit(wdata, paramsLS_diag{:}, 'replicates', 1);
                wgDs.params.sigma = wgDs.wv2sigma;
                
            case 'LS_diagw'
                tic
                % run ds wgmm
                [~, wgDsLs] = fitgmdist2wgmmLS(yds, gmmK, regVal, 3, dsgmmopt, dLow);
                % put diagonal along
                wgDsLsDiag = wgmm(wgDsLs.opts, wgDsLs.params);
                wgDsLsDiag.params.sigma = repmat(eye(prod(scAtlPatchSize)), [1,1,gmmK]); 
                wgDsLsDiag.params.W = wgDsLsDiag.params.sigma(:, 1:dLow, :);
                fprintf('ds GMM took %3.2f\n', toc);

                % run LS_diagw50 wgmm on top of ds wgmm.
                w = wtPatchesSmall(dataridx, :);
                w(w < wthr) = 0; 
%                 w = w ./ max(w(:)); % for algo
                
                tic
                wdata = struct('Y', yds, 'W', w, 'K', gmmK);
                paramsLS_diagw50 = [lsopts, 'init', 'wgmm', 'initArgs', struct('wgmm', wgDsLsDiag), 'modelName', 'wLatentSubspace'];
                wgDs = wgmmfit(wdata, paramsLS_diagw50{:}, 'replicates', 1, 'maxIter', 15);
                wgDs.params.sigma = wgDs.wv2sigma;
                fprintf('LS_WGMM took %3.2f\n', toc);
                
            otherwise
                error('unknwon wgmm method');
        end
        clear mxPatchesSmall wtPatchesSmall;
    
        % mark the previous wg
        wgDsSave{si} = wgDs;
        if si > 1
            wgDs.stats(1).wgprev = wgDsSave{si-1};
        end

        if si < numel(scSizes)
            tic
            % run patch reconstruction in atlas space
            % TODO: This is not great, since we'd like to avoid reconstructing with ds...
            wgDs.opts.verbose = 0;
            % reconstruct volumes
            subvolSize = size(vols.mxSubvolsSmall(:,:,:,1));
            reconVols = zeros(size(vols.mxSubvolsSmall));
            for vi = 1:size(vols.mxSubvolsSmall, 4)
                % reconstruct patches
                mxPatchesSmallVol = patchlib.vol2lib(vols.mxSubvolsSmall(:,:,:,vi), scAtlPatchSize);
                wtPatchesSmallVol = patchlib.vol2lib(vols.wtSubvolsSmall(:,:,:,vi), scAtlPatchSize);
                wdata = struct('Y', mxPatchesSmallVol, 'W', double(wtPatchesSmallVol > wthr), 'wts', wtPatchesSmallVol, 'K', gmmK);
                
                % quilt subvolumes
                blurMixDsPatchesSmallRecon = wgDs.recon(wdata, reconModel);
                reconVols(:,:,:,vi) = patchlib.quilt(blurMixDsPatchesSmallRecon, ...
                    subvolSize-scAtlPatchSize+1, scAtlPatchSize);
            end    
            wgDs.opts.verbose = 2;      
            fprintf('vol recon took %3.2f\n', toc);

            tic
            % prepare top-down volumes
            scAtlPatchSize = scSizes(si+1);
            blurSigma = 1/3 * patchSize ./ minPathSize;
            newvols = hierarchyPrepVolumesSimple(dsSubvols, wtSubvols, ...
                blurSigma, blurWindow, diffPad, patchSize, scAtlPatchSize);

            % resize reconstructed volumes up
            reconVolsUp = volresize(reconVols, size(newvols.mxSubvolsSmall));
            diffPadSubvol = size(newvols.mxSubvolsSmall) - size(newvols.mxSubvolSmallCrop);
            diffPadSubvol = diffPadSubvol(1:3)/2;
            reconVolsUp = cropVolume(reconVolsUp, [diffPadSubvol+1, 1], [size(newvols.mxSubvolsSmall(:,:,:,1)) - diffPadSubvol, nSubj]);

            % compute new volumes
            %w = min(newvols.wtSubvolSmallCrop .* wthrorig ./ prctile(newvols.wtSubvolSmallCrop(:), 100-100*percKnown), 1);
            wv = newvols.wtSubvolSmallCrop;
            mxSubvolSmallCrop = newvols.mxSubvolSmallCrop .* wv + reconVolsUp .* (1 - wv);

            % visualize
%             vidshow = 1;
%             view3Dopt(newvols.isoSubvolSmallCrop(:,:,:,vidshow), newvols.mxSubvolSmallCrop(:,:,:,vidshow), ...
%                 w(:,:,:,vidshow), reconVolsUp(:,:,:,vidshow), mxSubvolSmallCrop(:,:,:,vidshow));
%             drawnow;

            vols = newvols;
            vols.mxSubvolSmallCrop = mxSubvolSmallCrop;
            fprintf('vol extraction and merge took %3.2f\n', toc);
        end
    end
    
    if nargout > 1
        warning('NOTE: vols returned is from the previous (not final) iteration!.');
    end
    
    wg = wgDs;
    wg.params.sigma = wg.wv2sigma;
    save(wgmmMat, 'wg', '-v7.3'); 

end



function vols = hierarchyPrepVolumesSimple(dsSubvols, wtSubvols, ...
    blurSigma, blurWindow, diffPad, atlPatchSize, scAtlPatchSize)
% modelled after hierarchyPrepVolumes(), but does less stuff.

    % some processing
    blurSigma4D = [blurSigma, 0.001];
    blurWindow4D = [blurWindow, 1];
    volSize = size(dsSubvols);
    nSubj = volSize(4);
    volSize = volSize(1:3);
    ratio = scAtlPatchSize ./ atlPatchSize;

    % avoid nans
    wtSubvols = max(wtSubvols, 1e-7);

    if all(atlPatchSize == scAtlPatchSize) 
        % if max volume, don't blur, downsample, etc. just get patches.
        mxSubvolsSmall = dsSubvols;
        wtSubvolsSmall = wtSubvols;
        
        mxSubvolSmallCrop = dsSubvols;
        wtSubvolSmallCrop = wtSubvols;
    else
        % blur subvolumes
        blurDsWtSubvols = volblur(dsSubvols .* wtSubvols, blurSigma4D, blurWindow4D);
        blurWtSubvols = volblur(wtSubvols, blurSigma4D, blurWindow4D);

        % downsample
        subvolSmallSize = round(volSize .* ratio);
        mxSubvolsSmall = volresizeSimple(blurDsWtSubvols ./ blurWtSubvols, [subvolSmallSize, nSubj], 'simplelinear');
        wtSubvolsSmall = volresizeSimple(blurWtSubvols, [subvolSmallSize, nSubj], 'simplelinear');
        assert(isclean(mxSubvolsSmall));

        % crop
        diffPadSubvol = round(diffPad .* ratio);
        mxSubvolSmallCrop = cropVolume(mxSubvolsSmall, [diffPadSubvol+1, 1], [subvolSmallSize - diffPadSubvol, nSubj]);
        wtSubvolSmallCrop = cropVolume(wtSubvolsSmall, [diffPadSubvol+1, 1], [subvolSmallSize - diffPadSubvol, nSubj]);
    end
    
    vols = struct();
    vols.mxSubvolSmallCrop = mxSubvolSmallCrop;
    vols.wtSubvolSmallCrop = wtSubvolSmallCrop;
    vols.mxSubvolsSmall = mxSubvolsSmall;
    vols.wtSubvolsSmall = wtSubvolsSmall;
end
