function [logp, dsAtlPatchMeans] = subvolLogp(gmm, subvolLoc, subvolSize, atlPatchSize, dsSubjInAtlVol, dsSubjInAtlMaskVol)
% get the logp of all patches within a subvolume and return the quiltedSubvol
% along with the patch means
%
% Currently uses atlas-space heuristic to get the optimal cluster and the conditional rotation
% reconstruction method.

    % prepare libraries of patches of the subvolume and the mask-subvolume
    % TODO: note, we are extracting the mean from each patch, assuming that's how the gmms were
    % computed.
    dsAtlSubvol = cropVolume(dsSubjInAtlVol, subvolLoc, subvolLoc + subvolSize - 1);
    dsAtlPatches = patchlib.vol2lib(dsAtlSubvol, atlPatchSize);
    dsAtlPatchMeans = mean(dsAtlPatches, 2);
    dsAtlPatches = bsxfun(@minus, dsAtlPatches, dsAtlPatchMeans);
    
    dsAtlMaskSubvol = cropVolume(dsSubjInAtlMaskVol, subvolLoc, subvolLoc + subvolSize - 1);
    dsAtlMaskPatches = patchlib.vol2lib(dsAtlMaskSubvol, atlPatchSize);
    
    % TODO: just call wgmm.logp? it would include the prior
    
    % choose optimal cluster using the in-atlas heuristic measure
    K = size(gmm.mu, 1);
    logp = zeros(size(dsAtlPatches, 1), K);
    for k = 1:K
        atlMu = gmm.mu(k, :);
        atlSigma = gmm.sigma(:, :, k);
        
        atlmusel = bsxfun(@times, atlMu, dsAtlMaskPatches);
        atlpatchessel = dsAtlPatches .* dsAtlMaskPatches;
        logp(:, k) = logmvnpdf(atlpatchessel, atlmusel, atlSigma);
    end
end
