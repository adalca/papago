function [patches, layeridx, volidx] = matfile2patchcol(matfiles, vol, patchSize, loc, pad, isnii)
% MATFILE2PATCHCOL extract patch column from volumes in matfiles
%
% patches = matfile2patchcol(matfiles, vol, patchSize, loc, pad). See detailed discussion in
% vols2patchcol(). Here, matfiles is a (cell array of) matfile pointers, and vol is a string of the
% volume name in those files. This file is useful, compared to vols2patchcol for example, when there
% are too many volumes to hold in memory, or loading the volumes is quite slow. You still need, with
% this implementation, enough memory to hold all the subvolumes and the patch column.
%
% [patches, layeridx, volidx] = matfile2patchcol(...) also return layeridx and volidx (see
% subvols2patchcol()).
%
% See Also: vols2patchcol, subvols2patchcol, nii2patchcol, md2patchcol
%
% In the future, we would like to support this variant, but right now this is hard due to the way
% matfiles work:
% patches = matfile2patchcol(matfiles, vol, patchSize, loc, pad, isnii) allows for the specification
% on weather the volume vol is actually a nifti volume struct, in which case matfiles{i}.(vol).img
% will be used as the volume from which we extract the subvolume
%
%
%
% Contact: adalca@csail

    % input checking
    narginchk(5, 6);
    if ~iscell(matfiles)
        matfiles = {matfiles};
    end
    
    if ~exist('isnii', 'var')
        isnii = false;
    end
    if isnii
        warning('currently, loading a nii from inside the matfile does not use partial loading');
    end

    % extract subvolumes
    volRange = arrayfunc(@(x, d, p) (x - d: x + d + p - 1)', loc, pad, patchSize);
    subvols = cellfunc(@(ns) volcrop(ns, vol, volRange, isnii), matfiles);
    
    % get the patches 
    [patches, layeridx, volidx] = subspacetools.subvols2patchcol(subvols, patchSize);
end

function cvol = volcrop(matfile, vol, volRange, isnii)
    volRange = cellfunc(@(x) max(x, 1), volRange);
    
    volRange = cellfunc(@(x, y) min(x, y), volRange, mat2cellsplit(size(matfile, vol)));
    
    if ~isnii
        cvol = matfile.(vol)(volRange{:});
    else
        nii = matfile.(vol);
        cvol = nii.img(volRange{:});
    end
end
