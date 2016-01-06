function [patches, layeridx, volidx] = nii2patchcol(nii, patchSize, loc, pad)
% NII2PATCHCOL extract patch column from niftis
%
% patches = niis2patchcol(vols, patchSize, loc, pad). See detailed discussion in vols2patchcol().
% Here, niis is the (cell array of) inital, large nifti volumes from which we want to extract patch
% columns, or strings pointing to nifti files.
%
% [patches, layeridx, volidx] = niis2patchcol(...) also return layeridx and volidx (see
% vols2patchcol()).
%
% See Also: vols2patchcol, subvols2patchcol, md2patchcol, matfile2patchcol
%
% Contact: adalca@csail

    if ~iscell(nii)
        nii = {nii};
    end
    if ischar(nii{1}), nii = cellfunc(@loadNii, nii); end
    
    vols = cellfunc(@(ns) ns.img, nii);
    [patches, layeridx, volidx] = subspacetools.vols2patchcol(vols, patchSize, loc, pad);
end
