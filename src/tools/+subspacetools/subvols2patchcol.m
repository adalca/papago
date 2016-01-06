function [patches, layeridx, volidx] = subvols2patchcol(vols, patchSize)
% SUBVOLS2PATCHCOL extract patch column from subvolumes
%
% patches = subvols2patchcol(subvols, patchSize) extract the patches from the given volumes (cell of
% numerical arrays). in this setup, subvols2patchcol behaves like patchlib.vol2lib, except the
% patches are a vertical concatenation of the libraries you would get from vol2lib if vols is a
% cell.
%
% This function is meant to be called in the context of extracting a patch "column" our of several
% volumes. Specifically, in the simplest case, these subvolumes should be the patches themselves
% shaped in 3D -- in which case, all that this function does is reformat and stack the patches --
% hence a patch "column" However, in general, we expect the subvolumes will include the main patch,
% but also some patches around it. For example, if the patchSize is [9x9] but we want two extra
% patch in each direction each subvolume might be 13x13 (giving a total of 25 patches per volume).
% We assume the subvols were cropped to hold exactly the desired patches -- i.e. each subvolume only
% includes the middle of the desired patch column, plus the padding patches. 
%
% [patches, layeridx, volidx] = subvols2patchcol(...) also returns the layeridx and volidx. if
% patches is NxD, then layeridx and volidx are Nx1 each. Each entry of layeridx includes that
% patch's integer "layer" in terms of how far they are from the center patch. So the center patch
% will have layeridx 0, the next surrounding patches will have layeridx 1, etc. volidx indicates
% which volume each patch came from.
%
% TODO: see//erase libs2patchcells in volPatchImpute.
%
% See Also: vols2patchcol, nii2patchcol, md2patchcol, matfile2patchcol
%
% Contact: adalca@csail

    % check inputs
    narginchk(2, 2);
    if ~iscell(vols)
        vols = {vols};
    end
    % assert(all(isodd(patchSize)));
    vs = cellfunc(@(vol) size(vol) - patchSize, vols); vs = cat(1, vs{:});
    % assert(all(iseven(vs(:))));
    % TODO: check that all subvolumes are the same size
    
    % get libs
    [libs, locidx, ~, ~, volidx] = patchlib.vol2lib(vols, patchSize);
    patches = cat(1, libs{:});
    volidx = cat(1, volidx{:});
    
    % prepare layeridx
    cenloc = cellfunc(@(vol) ((size(vol) - patchSize + 1) + 1) / 2, vols);
    subs = cellfunc(@(l, v) ind2subvec(size(v), l), locidx, vols);
    layeridxc = cellfunc(@(x, loc) max(abs(bsxfun(@minus, x, loc)), [], 2), subs, cenloc);
    layeridx = cat(1, layeridxc{:});
end
