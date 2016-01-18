function [quiltedSubvol, minSubvolLoc, cntvol] = subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, ...
    crmethod, keepk, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin)
% reconstruct all patches within a subvolume and return the quiltedSubvol
% along with its location in the full volume (minSubvolLoc) and the number
% of patches that were used in reconstruting each pixel in the
% quiltedSubvol (cntvol)
%
% TODO: rename to subvolPatchRecon
% 
% subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepk, ...
%   dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubj, dsSubjWeight, atlLoc2SubjSpace, subjLoc2AtlSpace)
%
% subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepk, ...
%   dsSubjInAtlNii, dsSubjInAtlMaskNii, dsSubjNii, dsSubjWeightNii, tform)
%   
% Currently uses atlas-space heuristic to get the optimal cluster and the conditional rotation
% reconstruction method.

    % need: subjVol, subjWeightVol, crmethod, subjLoc2AtlSpace or regVal
    [dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg] ...
        = parseInputs(crmethod, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin{:});

    % prepare atlas locations to loop over
    atlLocs = patchlib.grid(subvolSize, atlPatchSize, 'sliding', 'sub');
    atlLocs = cellfunc(@(x) x(:), atlLocs);
    atlLocs = cat(2, atlLocs{:});
    
    % prepare reconstructions
    nLocs = size(atlLocs, 1);
    reconPatches = cell(nLocs, 1);
    reconLocs = cell(nLocs, 1);
    
    % go through each patch and reconstruct
    for i = 1:nLocs
        atlLoc = subvolLoc + atlLocs(i, :);

        % extract atlas patch
        dsAtlPatch = cropVolume(dsSubjInAtlVol, atlLoc, atlLoc + atlPatchSize - 1);
        dsAtlMaskPatch = cropVolume(dsSubjInAtlMaskVol, atlLoc, atlLoc + atlPatchSize - 1);
        
        % choose optimal cluster using the in-atlas heuristic measure
        K = size(gmm.mu, 1);
        logp = zeros(1, K);
        for k = 1:K
            atlMu = gmm.mu(k, :)';
            atlSigma = gmm.Sigma(:, :, k);
            logp(k) = logmvnpdf(dsAtlPatch(:)' .* dsAtlMaskPatch(:)', atlMu(:)' .* dsAtlMaskPatch(:)', atlSigma);
        end
        
        % get the optimal cluster via posteriors
        logpost = log(gmm.ComponentProportion) + logp;
        [~, optk] = sort(logpost, 'descend');
        
        % reconstruct in subject space
        for k = 1:keepk
            atlMu = gmm.mu(optk(k), :)';
            atlSigma = gmm.Sigma(:, :, optk(k));
            [reconPatches{i, k}, reconLocs{i, k}] = paffine.recon(atlMu, atlSigma, atlLoc, ...
                atlPatchSize, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArg);
        end
    end 

    
    % determine what region of the subvolume is quilted using the reconPatches
    reconLocsEnd = cellfunc(@(x, y) x  + size(y), reconLocs, reconPatches); 
    minSubvolLoc = min(reshape([reconLocs{:}], [3 length(reconLocs)]), [],2)';
    maxSubvolLoc = max(reshape([reconLocsEnd{:}], [3 length(reconLocs)]), [], 2)';
    
    % determine the size of the subvolume that will be quilted
    subvolSize = maxSubvolLoc - minSubvolLoc + 1; 

    % quilt the irregularly sized patches
    modReconLocs = cellfunc(@(x) x - minSubvolLoc + 1, reconLocs);
    [quiltedSubvol, cntvol] = patchlib.quiltIrregularPatches(modReconLocs, reconPatches, 'volSize', subvolSize);
    
end

function [subjAtlVol, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg] ...
    = parseInputs(crmethod, dsSubjInAtl, dsSubjInAtlMask, dsSubj, dsSubjWeight, varargin)

    if iscell(varargin{1})
        subjAtlVol = dsSubjInAtl;
        dsSubjInAtlMaskVol = dsSubjInAtlMask;
        dsSubjVol = dsSubj;
        dsSubjWeightVol = dsSubjWeight;
        atlLoc2SubjSpace = varargin{1};
        
        extraReconArg = [];
        if numel(varargin) > 1
            extraReconArg = varargin{2};
        end 
        
    else
        dsSubjInAtlNii = dsSubjInAtl;
        subjAtlVol = dsSubjInAtlNii.img;
        dsSubjInAtlMaskVol = dsSubjInAtlMask;
        dsSubjNii = dsSubj;
        dsSubjVol = dsSubjNii.img;
        dsSubjWeightVol = dsSubjWeight;
        
        tform = varargin{1};
        regVal = 0;
        if numel(varargin) > 1
            regVal = varargin{2};
        end 
        
        [subjLoc2AtlSpace, atlLoc2SubjSpace] = nii2cor3d(tform, dsSubjNii, dsSubjInAtlNii);
        extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);
    end
end




