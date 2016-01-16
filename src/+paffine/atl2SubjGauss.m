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
    
    % if need to compute R, compute it.
    needr = nargin > 5;
    if needr
        [R, subjInterpMask, subjPatchMins, subjPatchSize] = ...
            paffine.volCor2patchInterpmat(method, varargin{:});
    end
    

    % determine direction to execute
    switch method
        case {'forward', 'forwardModeling', 'forwardModel'}
            % get subject Gaussian parameters from atlas Gaussian Parameters
            subjSigma = pinv(R) * atlSigma * pinv(R');
            subjMu = pinv(R) * atlMu(:); 

            % regularize if desired
            if (needr && numel(varargin) == 4) || (~needr && numel(varargin) == 2)
                regVal = varargin{end};
                subjSigma = subjSigma + eye(size(subjSigma)) .* regVal;
            end

        case {'inverse', 'inverseModeling', 'inverseModel'}
            % go from atlas to subject via R
            subjSigma = R * atlSigma * R';
            subjMu = R * atlMu(:);
    
        otherwise
            error('not a valid method for sigma estimation');
    end
    
    % prepare extra outputs
    if nargout > 2
        assert(needr, 'can only ask for more outputs when computing R');
        varargout = {subjInterpMask, subjPatchMins, subjPatchSize};
        varargout = varargout(1:(nargout - 2));
    end
end
