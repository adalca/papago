function [patches, layeridx, volidx] = vols2patchcol(vols, patchSize, loc, pad)
% VOLS2PATCHCOL extract patch column from volumes
%
% patches = vols2patchcol(vols, patchSize, loc, pad). See detailed discussion in subvols2patchcol().
% Here, vols is the (cell array of) inital, large volumes from which we want to extract patch
% columns. patchSize is the size of the patch, loc is the center location of the main (center)
% patch. pad is the number of extra patches in each direction to include in the column. patchSize,
% loc and pad are the same size [1 x D]. if pad is [1, 1], for example, you'll get 9 patches per
% volume.
%
% [patches, layeridx, volidx] = vols2patchcol(...) also return layeridx and volidx (see
% subvols2patchcol()).
%
% See Also: subvols2patchcol, nii2patchcol, md2patchcol
%
% TODO: see//erase libs2patchcells in volPatchImpute.
%
% Contact: adalca@csail

    if ~iscell(vols)
        vols = {vols};
    end

    % extract subvolumes
    volRange = arrayfunc(@(x, d, p) (x - d: x + d + p - 1)', loc, pad, patchSize);
    subvols = cellfunc(@(ns) volcrop(ns, volRange), vols);
    
    % get the patches 
    [patches, layeridx, volidx] = subspacetools.subvols2patchcol(subvols, patchSize);
end

function cvol = volcrop(vol, volRange)
    volRange = cellfunc(@(x) max(x, 1), volRange);
    volRange = cellfunc(@(x, y) min(x, y), volRange, mat2cellsplit(size(vol)));
    cvol = vol(volRange{:});
end
