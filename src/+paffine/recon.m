function [reconPatch, logp, subjPatchMins] = recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
    subjVol, subjWeightVol, atlLoc2SubjSpace, method, varargin)
% reconstruct a subject patch given the atlas
%
% varargin is:
%   <regVal> if method is forward 
%   subjLoc2AtlSpace if method is inverse
%
% srcLoc2TgtSpace is computed via
%   srcLoc2TgtSpace = tform2cor3d(subj2Atl, size(subjVol), srcVoxDims, tgtVolSize, tgtVoxDims, dirn);
%
% TODO: use cropVolume to extract patch.
% TODO: allow just tform as opposed to cell map?

    % get subject gaussian coordinates and information
    [subjMu, subjSigma, subjInterpMask, subjPatchMins, subjPatchSize] = ...
        paffine.atl2SubjGauss(atlMu, atlSigma, method, atlLoc, atlPatchSize, atlLoc2SubjSpace, varargin{:});
           
    % reconstruct the patch
    subjPatch = cropVolume(subjVol, subjPatchMins, subjPatchMins + subjPatchSize - 1);
    subjWeightPatch = cropVolume(subjWeightVol, subjPatchMins, subjPatchMins + subjPatchSize - 1);
    [reconPatch, invBb] = paffine.reconSubjPatch(subjPatch, subjWeightPatch, subjInterpMask, subjMu, subjSigma);
    
    % compute the logp
    if nargout > 1
        logp = paffine.logpSubjPatch(subjPatch, subjWeightPatch, subjInterpMask, subjMu, subjSigma, invBb);
    end
end
