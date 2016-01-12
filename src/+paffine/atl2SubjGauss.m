function [subjMu, subjSigma] = atl2SubjGauss(atlMu, atlSigma, R, method, regVal)
% transform gaussian parameters from atlas space to subject space using a interpolation matrix R
% computed via cor2interpmat
%
% atlMu is a N-by-1
% atlSigma is a N-by-N matrix
% R is M-by-N if inverse (taking atlas to subject space), or N-by-M if forward (taking subject to
%   altas space)
% method is 'inverse' of 'forward'. 
% <regVal> is a diagonal regularization value for sigma: scalar of N-by-1.
%
% returnes subjMu and subjAtls are M-by-1 and M-by-M, respectively.
%
% part of papago


    narginchk(4, 5);
    
    switch method
        case {'forward', 'modelingFwd', 'forward-model'}
            subjSigma = pinv(R) * atlSigma * pinv(R');
            subjMu = pinv(R) * atlMu(:); 
            
            % regularize if desired
            if nargin >= 5
                subjSigma = subjSigma + eye(size(subjSigma)) .* regVal;
            end
            % newSigma = nearestSPD(subjSigma);
            
        case {'inverse', 'modelingInv', 'inverse-model'}
            subjSigma = R * atlSigma * R';
            subjMu = R * atlMu(:);
            
        otherwise
            error('not a valid method for sigma estimation');
    end