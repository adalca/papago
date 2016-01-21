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
    reconPatches = cell(nLocs, 1);
    reconLocs = cell(nLocs, 1);
    
    % extract relevant patches
    gatlLocs = bsxfun(@plus, subvolLoc, atlLocs);
    meanAtlPatch = zeros(nLocs, 1);
    dsAtlPatch = zeros(nLocs, prod(atlPatchSize));
    dsAtlMaskPatch = zeros(nLocs, prod(atlPatchSize));
    for i = 1:nLocs
        atlLoc = gatlLocs(i, :);

        % extract atlas patch and subtract mean
        dsAtlPatchtmp = cropVolume(dsSubjInAtlVol, atlLoc, atlLoc + atlPatchSize - 1);
        meanAtlPatch(i) = mean(dsAtlPatchtmp(:));
        dsAtlPatch(i, :) = dsAtlPatchtmp(:) - meanAtlPatch(i);
        dsAtlMaskPatchtmp = cropVolume(dsSubjInAtlMaskVol, atlLoc, atlLoc + atlPatchSize - 1);
        dsAtlMaskPatch(i, :) = dsAtlMaskPatchtmp(:);
    end
        
    % choose optimal cluster using the in-atlas heuristic measure
    K = size(gmm.mu, 1);
    logp = zeros(nLocs, K);
    for k = 1:K
        atlMu = gmm.mu(k, :)';
        atlSigma = gmm.sigma(:, :, k);
        atlmusel = bsxfun(@times, atlMu(:)', dsAtlMaskPatch);
        logp(:, k) = logmvnpdf(dsAtlPatch .* dsAtlMaskPatch, atlmusel, atlSigma);
    end
    
    for i = 1:nLocs
        % get the optimal cluster via posteriors
        logpost = log(gmm.pi) + logp(i, :);
        [~, optk] = sort(logpost, 'descend');
        
        % reconstruct in subject space
        for k = 1:keepk
            atlMu = gmm.mu(optk(k), :)' + meanAtlPatch(i);
            atlSigma = gmm.sigma(:, :, optk(k));
            [reconPatches{i, k}, reconLocs{i, k}] = paffine.recon(atlMu, atlSigma, gatlLocs(i, :), ...
                atlPatchSize, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArgs{:});
        end
        % reconPatches(i,:) = cellfunc(@(x) x + meanAtlPatch, reconPatches(i,:));    
    end

    
    % determine what region of the subvolume is quilted using the reconPatches
    reconLocsEnd = cellfunc(@(x, y) x  + size(y), reconLocs, reconPatches); 
    % minSubvolLoc = min(reshape([reconLocs{:}], [3 length(reconLocs)]), [],2)';
    minSubvolLoc = min(cat(1, reconLocs{:}), [], 1);
    maxSubvolLoc = max(cat(1, reconLocsEnd{:}), [], 1);
    % maxSubvolLoc = max(reshape([reconLocsEnd{:}], [3 length(reconLocs)]), [], 2)';
    
    % determine the size of the subvolume that will be quilted
    subvolSize = maxSubvolLoc - minSubvolLoc + 1; 

    % quilt the irregularly sized patches
    modReconLocs = cellfunc(@(x) x - minSubvolLoc + 1, reconLocs);
    [quiltedSubvol, cntvol] = patchlib.quiltIrregularPatches(modReconLocs, reconPatches, 'volSize', subvolSize);
    
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
