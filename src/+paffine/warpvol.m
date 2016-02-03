function regNii = warpvol(dsSubjNii, dsusSubjmasknii, interpSubjFile, atlSize, regOutFile)
% Transform sparse-slice volumes according to sparse-slice interpolant matrix, 
% as opposed to warping dsXusX images
%
% Naming: S for Subject Space; A for Atlas space; Ss for sparse-slice (small in z) subject space
%
% Steps:
% 1. Compute R that takes S to A, so R is |A|-by-|S|
% 2. Compute T that takes Ss to A, so T is |A|-by-|Ss| (a lot skinnier)
% 3. Transform subjects via T and save.
%
% Note, the loading of R replaces the following code, which requires the atl2subjTform
%   (subj2atlTform.invert)
% atl2subjTform = tform.invert;
% subjVolSize = size(dsusmasknii.img);
% subjVoxDims = dsusmasknii.hdr.dime.pixdim(2:4);
% atlVolSize = size(atlnii.img);
% atlVoxDims = atlnii.hdr.dime.pixdim(2:4);
% atl2subjCor = tform2cor3d(atl2subjTform, atlVolSize, atlVoxDims, subjVolSize, subjVoxDims);
% subj2atlR = cor2interpmat(subjVolSize, atl2subjCor);

    if ischar(dsSubjNii), dsSubjNii = loadNii(dsSubjNii); end
    if ischar(dsusSubjmasknii), dsusSubjmasknii = loadNii(dsusSubjmasknii); end
    if ischar(atlSize)
        if sys.isfile(atlSize)
            atlNii = loadNii(atlSize);
            atlSize = size(atlNii.img);
        else
            atlSize = str2num(atlSize);
        end
    end
    
    % load subj2atlR
    load(interpSubjFile, 'subj2atlR');  
    
    % extract the subset of R for the voxels we trust in the ds data
    mask = logical(dsusSubjmasknii.img);
    T = subj2atlR(:, mask(:));

    % prepare re-normalization factor given that we deleted some voxels
    normfactT = 1./sum(T, 2);

    % extract the downsampled subject volume
    dsvol = dsSubjNii.img(mask(:));

    % compute warped volume
    warpedVol = reshape(normfactT .* (T * dsvol(:)), atlSize); 
    
    % save modality
    if exist('regoutfile', 'var')
        regNii = dsusSubjmasknii;
        regNii.img = warpedVol;
        regNii.hdr.dime.dim(2:4) = atlSize;
        saveNii(regNii, regOutFile)
    end
end
