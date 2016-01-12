function [reconPatch, logp] = recon(atlMu, atlSigma, atlLoc, atlPatchSize, subjVol, subjWeightVol, ...
    atlLoc2SubjSpace, method, varargin)
% reconstruct a subject patch given the atlas
%
% srcLoc2TgtSpace is computed via
%   srcLoc2TgtSpace = tform2cor3d(subj2Atl, size(subjVol), srcVoxDims, tgtVolSize, tgtVoxDims, dirn);


    % get the subject space coordinates
    [subjPatchRange, subjPatchMins, subjPatchSize] = ...
        paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);

    switch method
        case {'forward', 'modelingFwd', 'forward-model'}
            regVal = varargin{1};
            
            % obtain R subj --> atl. R is |atl|-by-|subj|
            atlPatchRange = arrayfunc(@(a, p) a:p, atlLoc, atlPatchSize);
            atlLoc2SubjSpacePatch = extractAndNormalizePatchCor(atlLoc2SubjSpace, atlPatchRange, atlLoc);
            [subj2AtlR, subjMask, ~] = cor2interpmat(subjPatchSize, atlLoc2SubjSpacePatch);
            
            % get the subject-space gaussian parameters
            [subjMu, subjSigma] = paffine.atl2SubjGauss(atlMu, atlSigma, subj2AtlR, method, regVal);
            
        case {'inverse', 'modelingInv', 'inverse-model'}
            subjLoc2AtlSpace = varargin{1};
            
            % obtain R atl --> subj. R is |subj|-by-|atl|
            subjLoc2AtlSpacePatch = extractAndNormalizePatchCor(subjLoc2AtlSpace, subjPatchRange, subjPatchMins);
            [atl2SubjR, ~, subjMask] = cor2interpmat(atlPatchSize, subjLoc2AtlSpacePatch);
            
            % get the subject-space gaussian parameters
            [subjMu, subjSigma] = paffine.atl2SubjGauss(atlMu, atlSigma, atl2SubjR, method);
            
        otherwise
            error('paffine: unknown recon method');
    end
           
    % reconstruct the patch
    subjPatch = subjVol(subjPatchRange{:});
    subjWeightPatch = subjWeightVol(subjPatchRange{:});
    [reconPatch, invBb] = paffine.reconSubjPatch(subjPatch, subjWeightPatch, subjMask, subjMu, subjSigma);
    
    % compute the logp
    if nargout > 1
        logp = logpSubjPatch(subjPatch, subjWeightPatch, subjMask, subjMu, subjSigma, invBb);
    end
end

function patchCor = extractAndNormalizePatchCor(volCor, range, mins)
    patchCor = cellfunc(@(x) x(range), volCor);
    patchCor = cellfunc(@(x, m) x - m + 1, patchCor, mat2cellsplit(mins));
end
