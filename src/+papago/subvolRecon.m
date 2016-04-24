function [quiltedSubvol, minSubvolLoc, cntvol, modReconLocs, reconPatches] = ...
    subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepk, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin)
% reconstruct all patches within a subvolume and return the quiltedSubvol
% along with its location in the full volume (minSubvolLoc) and the number
% of patches that were used in reconstruting each pixel in the
% quiltedSubvol (cntvol)
%
% TODO: rename to subvolPatchRecon
% 
% subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepk, ...
%   dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubj, dsSubjWeight, atlLoc2SubjSpace, subjLoc2AtlSpace, <bigR>)
%
% subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepk, ...
%   dsSubjInAtlNii, dsSubjInAtlMaskNii, dsSubjNii, dsSubjWeightNii, tform, <bigR>)
%   
% Currently uses atlas-space heuristic to get the optimal cluster and the conditional rotation
% reconstruction method.

    % need: subjVol, subjWeightVol, crmethod, subjLoc2AtlSpace or regVal
    [dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArgs] ...
        = parseInputs(crmethod, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin{:});

    % prepare atlas locations to loop over
    atlLocs = patchlib.grid(subvolSize, atlPatchSize, 'sliding', 'sub');
    atlLocs = cellfunc(@(x) x(:), atlLocs);
    atlLocs = cat(2, atlLocs{:});
    
    % prepare reconstructions
    nLocs = size(atlLocs, 1);
    gatlLocs = bsxfun(@plus, subvolLoc - 1, atlLocs);
    reconPatches = cell(nLocs, 1);
    reconLocs = cell(nLocs, 1);
    reconWeights = cell(nLocs, 1);
    
    % logp computation for each patch, for each cluster
    [logp, meanAtlPatch] = papago.subvolLogp(gmm, subvolLoc, subvolSize, atlPatchSize, dsSubjInAtlVol, dsSubjInAtlMaskVol);
    
    % reconstruction
    for i = 1:nLocs
        % get the optimal cluster via posteriors
        logpost = log(gmm.pi) + logp(i, :);
        [~, optk] = sort(logpost, 'descend');
        explogpost = exp(logpost(optk) - max(logpost(optk)));
        
        % reconstruct in subject space
        for k = 1:keepk
            atlMu = gmm.mu(optk(k), :)' + meanAtlPatch(i);
            atlSigma = gmm.sigma(:, :, optk(k));
            [reconPatches{i, k}, reconLocs{i, k}, newSigmas{i,k}] = paffine.recon(atlMu, atlSigma, gatlLocs(i, :), ...
                atlPatchSize, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArgs{:});
            
            reconWeights{i, k} = reconPatches{i, k} * 0 + explogpost(optk(k));
%             reconWeights{i, k} = reconPatches{i, k} * 0 + logpost(optk(k));
        end
        % reconPatches(i,:) = cellfunc(@(x) x + meanAtlPatch, reconPatches(i,:));
    end

    % determine what region of the subvolume is quilted using the reconPatches
    reconLocsEnd = cellfunc(@(x, y) x  + size(y), reconLocs, reconPatches); 
    minSubvolLoc = min(cat(1, reconLocs{:}), [], 1);
    maxSubvolLoc = max(cat(1, reconLocsEnd{:}), [], 1);
    
    % determine the size of the subvolume that will be quilted
    subvolSize = maxSubvolLoc - minSubvolLoc + 1; 

    % quilt the irregularly sized patches
    modReconLocs = cellfunc(@(x) x - minSubvolLoc + 1, reconLocs);
%     [quiltedSubvol, cntvol] = patchlib.quiltIrregularPatches(modReconLocs, reconPatches, ...
%         'weightPatches', reconWeights, 'volSize', subvolSize);
    [quiltedSubvol, cntvol] = patchlib.quiltIrregularPatches(modReconLocs, reconPatches, ...
            'volSize', subvolSize);    
end

function [subjAtlVol, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArgs] ...
    = parseInputs(crmethod, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin)

    if iscell(varargin{1})
        subjAtlVol = dsSubjInAtl;
        dsSubjInAtlMaskVol = dsSubjInAtlMask;
        dsSubjVol = dsSubj;
        dsSubjWeightVol = dsSubjWeight;
        atlLoc2SubjSpace = varargin{1};
        
        extraReconArgs = {};
        if numel(varargin) > 1
            extraReconArgs = varargin(2:end);
        end 
        
    else
        dsSubjInAtlNii = dsSubjInAtl;
        subjAtlVol = dsSubjInAtlNii.img;
        dsSubjInAtlMaskVol = dsSubjInAtlMask;
        dsSubjNii = dsSubj;
        dsSubjVol = dsSubjNii.img;
        dsSubjWeightVol = dsSubjWeight;
        
        tform = varargin{1};
        
        [subjLoc2AtlSpace, atlLoc2SubjSpace] = nii2cor3d(tform, dsSubjNii, dsSubjInAtlNii);
        extraReconArgs = ifelse(strcmp(crmethod, 'inverse'), {subjLoc2AtlSpace, varargin{2:end}}, varargin(2:end));
    end
end
