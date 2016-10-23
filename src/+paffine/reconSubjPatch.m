function [reconPatch, invBb, varargout] = reconSubjPatch(subjPatch, subjWeightPatch, subjPatchValidRegion, subjMu, subjSigma)
% RECONSUBJPATCH reconstruct a subject-space patch
%
% subjPatch is the whole patch in subject space
% subjWeight is 0-1 weight patch in subject space
% subjPatchValidRegion is the 0-1 mask from the transform specifying the atlas patch shape inside
%   the subject patch. Think of this as the "shape" of the atlas patch rotated in subject space.
% subjMu is the nKnown-by-1 subject mean
% subjSigma is the nKnown-by-nKnown subject covariance
%
% reconPatch is the same size as subjPatch, with all the voxels inside the subjPatchValidRegion
% invBb is a potentially-useful matrix of inv(B) * b, where b is the vector of mean-centered known
% and valid voxels, and B is the covariance of those voxels.
    
    varargout = {};
    if nargout > 2
        varargout = {[]};
    end
    assert(islogical(subjPatchValidRegion));

    % get the part of the subject patch corresponding to the atlas patch
    subjPatchWithNans = subjPatch;
    subjPatchWithNans(~subjWeightPatch) = nan;
    subjPatchAtlVoxels = subjPatchWithNans(subjPatchValidRegion);
    
    % impute the voxels
    selSubjMu = subjMu(subjPatchValidRegion);
    selSubjSigma = subjSigma(subjPatchValidRegion, subjPatchValidRegion);
    [subjPatchAtlReconVoxels, invBb, varargout{:}] = ...
        inpaintWithGaussConditional(subjPatchAtlVoxels(:), selSubjMu, selSubjSigma);
    
    % recon the subject patch
    reconPatch = nan(size(subjPatch));
    reconPatch(subjPatchValidRegion) = subjPatchAtlReconVoxels;
