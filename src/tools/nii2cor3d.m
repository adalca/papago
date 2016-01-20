function [subjLoc2AtlSpace, varargout] = nii2cor3d(tform, dsSubjNii, dsSubjInAtlNii)
% see tform2cor3d

    if ischar(tform), load(tform, 'tform'); end
    if ischar(dsSubjNii), dsSubjNii = loadNii(dsSubjNii, tempname, false, true); end
    if ischar(dsSubjInAtlNii), dsSubjInAtlNii = loadNii(dsSubjInAtlNii, tempname, false, true); end
    
    % prepare necessary inputs for conditional-based reconstruction
    subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
    subjVolSize = dsSubjNii.hdr.dime.dim(2:4);
    atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
    atlVolSize = dsSubjInAtlNii.hdr.dime.dim(2:4);
    
    subjLoc2AtlSpace = tform2cor3d(tform, subjVolSize, subjDims, atlVolSize, atlDims);
    
    if nargout > 1 % atlLoc2SubjSpace
        varargout{1} = tform2cor3d(tform, subjVolSize, subjDims, atlVolSize, atlDims, 'backward');
    end
end
