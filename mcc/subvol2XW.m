function [X, W] = subvol2XW(subvolfile, volnames, subtractpatchmean, patchSize, nSamples, randomseed)


    % check inputs.
    narginchk(5, 6)
    if ischar(volnames), volnames = str2cell([volnames ','], ','); end
    if ischar(nSamples), nSamples = str2num(nSamples); end 
    if ischar(patchSize), patchSize = str2num(patchSize); end
    if ischar(subtractpatchmean), subtractpatchmean = str2num(subtractpatchmean); end
    if nargin >= 7 && ischar(randomseed), randomseed = str2num(randomseed); end
    if nargin < 7, randomseed = now; end

    % get subvolumes
    volumes = struct2cell(load(subvolfile, volnames{:}));
    volumes = dimsplit(2, volumes(:)');
    volSizes = cell2mat(cellfunc(@size, volumes{1}(:)));
    effVolSizes = bsxfun(@minus, volSizes, patchSize + 1); % effective (samplable) size
    
    % prepare sampling amounts. 
    if isscalar(nSamples) && nSamples == 0 % all
        nSamples = sum(prod(effVolSizes, 2));
    elseif isscalar(nSamples) && nSamples < 1 % fraction of all
        nSamples = sum(prod(effVolSizes, 2)) * nSamples;
    end
    
    % sample patches
    rng(randomseed); % set random seed for sampling.
    patches = patchlib.vol2samples(nSamples, patchSize, volumes{:});
    
    % prepare data (X,W) 
    if size(volumes, 2) == 2
        X = patches{1};
        W = max(patches{2}, params.minWeight);
    else
        X = patches;
        W = ones(size(X));
    end
    
    % subtract patch means if necessary
    if subtractpatchmean;
        X = bsxfun(@minus, X, mean(X, 2));
    end
    