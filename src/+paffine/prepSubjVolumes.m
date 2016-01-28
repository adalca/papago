function [atl2subjR , subj2atlR, subjLoc2AtlSpace, atlLoc2SubjSpace] = ...
    prepSubjVolumes(subjInAtlTform, dsSubjNii, dsSubjInAtlNii, outFilename)
% prepare correspondance and interpolation volumes for a given subject useful for the paffine
% transformations.
%
% all inputs can be filenames to the right file type

    % get the correspondances from nifti files.
    [subjLoc2AtlSpace, atlLoc2SubjSpace] = nii2cor3d(subjInAtlTform, dsSubjNii, dsSubjInAtlNii); 
    
    % prepare the interpolation matrices.
    [atl2subjR, ~, ~] = cor2interpmat(size(atlLoc2SubjSpace{1}), subjLoc2AtlSpace); 
    [subj2atlR, ~, ~] = cor2interpmat(size(subjLoc2AtlSpace{1}), atlLoc2SubjSpace);

    % save 
    save(outFilename, 'atl2subjR', 'subj2atlR', 'subjLoc2AtlSpace', 'atlLoc2SubjSpace'); 
    