function [subjMu, subjSigma, varargout] = atl2SubjGauss(atlMu, atlSigma, method, varargin)
% transform gaussian parameters from atlas space to subject space using a interpolation matrix R
% computed via cor2interpmat
%
% [subjMu, subjSigma] = atl2SubjGauss(atlMu, atlSigma, method, R, <regVal>)
%
% [subjMu, subjSigma] = atl2SubjGauss(atlMu, atlSigma, method, atlLoc, atlPatchSize, atlLoc2SubjSpace, regVal or subjLoc2AtlSpace)
%
% atlMu is a N-by-1
% atlSigma is a N-by-N matrix
% R is M-by-N if inverse (taking atlas to subject space), or N-by-M if forward (taking subject to
%   altas space)
% method is 'forward', 'forwardModeling', 'forwardModel', 'inverse', 'inverseModeling', 'inverseModel'
% <regVal> is a diagonal regularization value for sigma: scalar of N-by-1.
%
% returnes subjMu and subjAtls are M-by-1 and M-by-M, respectively.
%
% part of papago

    % input checking and output preparation
    narginchk(4, 7);
    
    varargout = cell(1, nargout - 2);

    % determine direction to execute
    switch method
        case {'forward', 'forwardModeling', 'forwardModel'}
            atl2SubjFcn = @atl2SubjFwd;

        case {'inverse', 'inverseModeling', 'inverseModel'}
            atl2SubjFcn = @atl2SubjBwd;
    
        otherwise
            error('not a valid method for sigma estimation');
    end
    
    % call transforming function
    [subjMu, subjSigma, varargout{:}] = atl2SubjFcn(atlMu, atlSigma, varargin{:});
end
    
function [subjMu, subjSigma, varargout] = atl2SubjFwd(atlMu, atlSigma, varargin)
% [subjMu, subjSigma] = atl2SubjFwd(atlMu, atlSigma, R, <regVal>)
% or
% [subjMu, subjSigma, subbjmask, subjPatchMins, subjPatchSize] = 
%   atl2SubjFwd(atlMu, atlSigma, atlLoc, atlPatchSize, atlLoc2SubjSpace, <regVal>)

    needr = nargin > 4;
    
    % obtain subj2AtlR: R subj --> atl. R is |atl|-by-|subj|
    if needr 
        warning('forward in currently untested!');
        [atlLoc, atlPatchSize, atlLoc2SubjSpace] = varargin{1:3};
        [~, subjPatchMins, subjPatchSize] = paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);

        atlPatchRange = arrayfunc(@(a, p) a:p, atlLoc, atlLoc + atlPatchSize - 1);
        atlLoc2SubjSpacePatch = extractAndNormalizePatchCor(atlLoc2SubjSpace, atlPatchRange, subjPatchMins);
        [R, subjMask, ~] = cor2interpmat(subjPatchSize, atlLoc2SubjSpacePatch);
        
    else
        R = varargin{:};
    end

    % get subject Gaussian parameters from atlas Gaussian Parameters
    subjSigma = pinv(R) * atlSigma * pinv(R');
    subjMu = pinv(R) * atlMu(:); 

    % regularize if desired
    if (needvar && numel(varargin) == 4) || (~needvar && numel(varargin) == 2)
        regVal = varargin{end};
        subjSigma = subjSigma + eye(size(subjSigma)) .* regVal;
    end
    % newSigma = nearestSPD(subjSigma);
    
    % prepare extra outputs
    if nargout > 2
        assert(needr, 'can only ask for more outputs when computing R');
        varargout = {subjMask, subjPatchMins, subjPatchSize};
        varargout = varargout(1:(nargout - 2));
    end
end

function [subjMu, subjSigma, varargout] = atl2SubjBwd(atlMu, atlSigma, varargin)
% [subjMu, subjSigma] = atl2SubjBwd(atlMu, atlSigma, R)
% or
% [subjMu, subjSigma, subbjmask, subjPatchMins, subjPatchSize] = 
%   atl2SubjBwd(atlMu, atlSigma, atlLoc, atlPatchSize, atlLoc2SubjSpace, subjLoc2AtlSpace)

    needr = nargin > 3;

    % obtain atl2SubjR: R atl --> subj. R is |subj|-by-|atl|
    if needr 
        [atlLoc, atlPatchSize, atlLoc2SubjSpace, subjLoc2AtlSpace] = varargin{:};
        [subjPatchRange, subjPatchMins, subjPatchSize] = paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);
        subjLoc2AtlSpacePatch = extractAndNormalizePatchCor(subjLoc2AtlSpace, subjPatchRange, atlLoc);
        [R, ~, subjMask] = cor2interpmat(atlPatchSize, subjLoc2AtlSpacePatch);
    else
        R = varargin{:};
    end

    % go from atlas to subject via R
    subjSigma = R * atlSigma * R';
    subjMu = R * atlMu(:);
    
    % prepare extra outputs
    if nargout > 2
        assert(needr, 'can only ask for more outputs when computing R');
        varargout = {subjMask, subjPatchMins, subjPatchSize};
        varargout = varargout(1:(nargout - 2));
    end
end
    
function patchCor = extractAndNormalizePatchCor(volCor, range, mins)
    patchCor = cellfunc(@(x) x(range{:}), volCor);
    patchCor = cellfunc(@(x, m) x - m + 1, patchCor, mat2cellsplit(mins));
end
