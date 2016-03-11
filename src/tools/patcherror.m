function err = patcherror(testpatches, truepatches, fun, preprocessmethod, varargin)
% PATCHERROR computer the patch error according to some error function. 
%   err = patcherror(testpatches, truepatches, fun) compares the stacked
%   test and true patches, each a N x prod(patchSize), according to the
%   comparison function fun. If not provided or if empty, fun is @msd.
%
%   err = patcherror(testpatches, truepatches, fun, function) allows for
%   some pre-processing of the given patches. Options are:
%       'flat' (default) - do no pre-processing
%       'subpatch'  - apply fun on subpatch. Need extra args: patchSize, subpatchstart, subpatchend
%       'gauss'     - msd on full patch weighted by a gaussian (patches are multiplied by gaussian blob)
%                   Need args: patchSize, sigma

    if nargin == 2 || isempty(fun)
        fun = @msd;
    end

    vst = @(x) x(:)';

    % prepare a processing method of the patch variables
    switch preprocessmethod
        case 'flat'
            processfun = @(x) x; % do nothing
            
        case 'gauss'
            gblob = gaussianblob(varargin{1}, varargin{2});
            processfun = @(x) (bsxfun(@times, x, gblob(:)')); % vst(...)
            
        case 'subpatch'
            patchSize = varargin{1};
            subpatchstart = varargin{2};
            subpatchend = varargin{3};
    
            % prepare extraction 
            extract1fun = @(x) vst(cropVolume(reshape(x, patchSize), subpatchstart, subpatchend));
            extractcellfun = @(x) cellfunc(extract1fun, dimsplit(1, x));
            catcell = @(n, c) cat(n, c{:});
            processfun = @(x) vst(catcell(1, extractcellfun(x)));
            
        otherwise
            % assume do nothing (flat)
            processfun = @(x) x; % do nothing
    end
    
    % compute the error
    err = fun(processfun(testpatches), processfun(truepatches));
end

function h = gaussianblob(patchSize, sigma) 
% gaussian blob of size patchSize

    range = arrayfunc(@(x) -x:x, (patchSize-1)/2);
    grid = cellfunc(@(x) x(:), ndgrid2cell(range{:}));
    grid = cat(2, grid{:});
    h = exp(-sum(grid.*grid, 2)/(2*sigma.^2));
end
