function [atlVolSize, subjmasks, dsregmaskvols, locVolumeAtlas, locVolumeSubjects, dsregvol] =  ...
    preloadTesting(reconSubjs, mods, testmdfile)
% PRELOADTESTING preload (in this case, actually load volumes) for testing.
%
% TODO: have some of the same options as before: load volumes or just "prepare" them ?
%
% TODO: should the loading actually pass in the md? probably...
%
% mods: 
%   dsmask sprintf('brainDs%dUs%dMask', dsRate, usRate)
%   ds sprintf('brainDs%dUs%d', dsRate, usRate)
%   dsreg sprintf('brainDs%dUs%dReg', dsRate, usRate)
%   dsregmask sprintf('brainDs%dUs%dRegMask', dsRate, usRate)
%   dsregmat sprintf('brainDs%dUs%dRegMat', dsRate, dsRate)


    testmd = load(testmdfile);

    
    for s = 1:numel(reconSubjs)
        reconSubj = reconSubjs{s};

        % load subject data in subject space
        % necessary for gaussian conditional reconstruction
        subjdsnii = testmd.md.loadModality(mods.ds, reconSubj);
        subjdsvol = double(subjdsnii.img);
        subjVolSize = size(subjdsvol);
        subjDims = subjdsnii.hdr.dime.pixdim(2:4);
        rSubjSpace = imref3d(subjVolSize, subjDims(2), subjDims(1), subjDims(3));
        subjmasks{s} = nii2vol(testmd.md.loadModality(mods.dsmask, reconSubj));

        % load subject data in atlas space (i.e. registered data)  
        % necessary for gaussian conditional reconstruction
        dsregnii = testmd.md.loadModality(mods.dsreg, reconSubj);
        dsregvol = dsregnii.img;
        atlVolSize = size(dsregvol); % size of volumes in atlas space
        atlDims = dsregnii.hdr.dime.pixdim(2:4);
        rAtlSpace = imref3d(atlVolSize, atlDims(2), atlDims(1), atlDims(3));
        dsmaskregnii = testmd.md.loadModality(mods.dsregmask, reconSubj);
        dsregmaskvols{s} = dsmaskregnii.img;
        regtformmat = load(testmd.md.getModality(mods.dsregmat, reconSubj));

        % get the locations in subject volume that correspond to voxel locations in the atlas volume 
        locVolumeAtlas{s} = getCorrespondingLoc(atlVolSize, regtformmat.tform, rAtlSpace, rSubjSpace, subjVolSize);
        locVolumeSubjects{s} = getCorrespondingLoc(subjVolSize, regtformmat.tform.invert, rSubjSpace, rAtlSpace, atlVolSize);
    end
    