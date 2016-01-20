function [atl2subjR , subjLoc2AtlSpace, atlLoc2SubjSpace] = prepSubjHelpfulVolumes(subjInAtlTform, dsSubjNii, dsSubjInAtlNii, outFilename)


[subjLoc2AtlSpace, atlLoc2SubjSpace] = nii2cor3d(subjInAtlTform, dsSubjNii, dsSubjInAtlNii); 
[atl2subjR, ~, ~] = cor2interpmat(size(atlLoc2SubjSpace{1}), subjLoc2AtlSpace); 

save(outFilename, 'atl2subjR', 'subjLoc2AtlSpace', 'atlLoc2SubjSpace'); 