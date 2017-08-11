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
    % dsAtlPatchMeans = mean(dsAtlPatches, 2);
    dsAtlPatchMeans = zeros(size(dsAtlPatches, 1), 1);
    warning('subtractMean forced FALSE!. Fix this elegantly!');
    dsAtlPatches = bsxfun(@minus, dsAtlPatches, dsAtlPatchMeans);
    
    dsAtlMaskSubvol = cropVolume(dsSubjInAtlMaskVol, subvolLoc, subvolLoc + subvolSize - 1);
    dsAtlMaskPatches = patchlib.vol2lib(dsAtlMaskSubvol, atlPatchSize);
    
    % TODO: just call wgmm.logp? it would include the prior
    
    % choose optimal cluster using the in-atlas heuristic measure
    K = size(gmm.params.mu, 1);
    logp = zeros(size(dsAtlPatches, 1), K);
    for k = 1:K
        atlMu = gmm.params.mu(k, :);
        if ~isfield(gmm.params, 'sigma') || size(gmm.params.sigma, 3) < k
            gmm.params.sigma(:, :, k) = gmm.params.W(:, :, k) * gmm.params.W(:, :, k)' + gmm.params.sigmasq(k) * eye(size(gmm.params.W, 1));
        end
            
        atlSigma = gmm.params.sigma(:, :, k);
            
        
        atlmusel = bsxfun(@times, atlMu, dsAtlMaskPatches);
        atlpatchessel = dsAtlPatches .* dsAtlMaskPatches;
        logp(:, k) = logmvnpdf(atlpatchessel, atlmusel, atlSigma);
    end
    
%     warning('New logp estimation!!')
%     logp = zeros(size(dsAtlPatches, 1), K);
%     for k = 1:K
%         atlMu = gmm.mu(k, :);
%         atlSigma = gmm.sigma(:, :, k);
%         
%         for i = 1:size(dsAtlPatches, 1)
%             idx = dsAtlMaskPatches(i, :) > 0.5;
%             atlmusel = atlMu(idx);
%             atlSigmasel = atlSigma(idx, idx);
%             atlpatchessel = dsAtlPatches(i, idx);
%             logp(i, k) = logmvnpdf(atlpatchessel, atlmusel, atlSigmasel);
%         end
%     end
end
