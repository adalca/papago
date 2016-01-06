function [subgridlocs, subGridLength, subGridSize] = prepgrid(atlVolSize, patchSize, patchOverlap, reconVolRange)

    if nargin == 3
        reconVolRange = {atlVolSize * 0, atlVolSize+1};
    end

    % obtain the grid points to operate on, given patchsize and overlap
    [gridsub, newVolSize, gridsize, overlap] = patchlib.grid(atlVolSize, patchSize, patchOverlap, 'sub');
    gridsubc = cellfunc(@(x) x(:), gridsub);
    gridlocs = cat(2, gridsubc{:});

    % compute which grid sub-points to run on by comparing gridsub and reconVolRange
    subvolgrididx = find(all(bsxfun(@gt, gridlocs, reconVolRange{1}) & bsxfun(@lt, gridlocs, reconVolRange{2}), 2));
    subgridlocs = gridlocs(subvolgrididx, :);
    subGridLength = numel(subvolgrididx);
    
    % reshape
    subGridSize = cellfun(@numel, cellfunc(@unique, dimsplit(2, subgridlocs)));
    assert(prod(subGridSize) == subGridLength);