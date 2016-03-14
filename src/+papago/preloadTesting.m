function [atlVolSize, subjmasks, dsregmaskvols, atlLoc2SubjSpace, subjLoc2AtlSpace, dsregvol] =  ...
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
        subjmasks{s} = nii2vol(testmd.md.loadModality(mods.dsmask, reconSubj));

        % load subject data in atlas space (i.e. registered data)  
        % necessary for gaussian conditional reconstruction
        dsregnii = testmd.md.loadModality(mods.dsreg, reconSubj);
        dsregvol = dsregnii.img;
        atlVolSize = size(dsregvol); % size of volumes in atlas space
        atlDims = dsregnii.hdr.dime.pixdim(2:4);
        dsmaskregnii = testmd.md.loadModality(mods.dsregmask, reconSubj);
        dsregmaskvols{s} = dsmaskregnii.img;
        regtformmat = load(testmd.md.getModality(mods.dsregmat, reconSubj));

        % get the locations in subject volume that correspond to voxel locations in the atlas volume 
        subjLoc2AtlSpace{s} = tform2cor3d(regtformmat.tform, subjVolSize, subjDims, atlVolSize, atlDims);
        atlLoc2SubjSpace{s} = tform2cor3d(regtformmat.tform, subjVolSize, subjDims, atlVolSize, atlDims, 'backward');
    end
    