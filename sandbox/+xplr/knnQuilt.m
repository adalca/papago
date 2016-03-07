
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = [52 43 34]; % LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
reconSubj = 1; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'adni'; warning('I dont have working bukner volumes for some reason'); 

% recon params
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

% modalities
ds = 5;
us = 5; warning('only works for us5 right now'); 


isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
dsSubjInAtlMod = sprintf('brainDs%dUs%dReg', ds, us);
dsSubjInAtlMaskMod = sprintf('brainDs%dUs%dRegMask', ds, us);
dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);

dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);

% isoSubjInAtlMod = sprintf('Iso2Ds%dUs%dsizeReg', ds, us);
% dsSubjInAtlMod = sprintf('Ds%dUs%dReg', ds, us);
% dsSubjInAtlMaskMod = sprintf('Ds%dUs%dRegMask', ds, us);
% dsSubjInAtlMatMod = sprintf('Ds%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
% dsInterpSubjInAtlMod = sprintf('Ds%dUs%dInterpReg', ds, us);
% 
% dsSubjMod = sprintf('Ds%dUs%d', ds, us);
% dsSubjMaskMod = sprintf('Ds%dUs%dMask', ds, us);
% isoSubjMod = sprintf('Iso2Ds%dUs%dsize', ds, us);


%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load various data columns from training set
fnames = fullfile(SYNTHESIS_DATA_PATH, traindataset, 'md', [sys.usrname, '_restor_md_*']);
trainmd = loadmd(fnames);
[bucknerIsoPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(trainmd, isoSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsMaskPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMaskMod, atlPatchSize, atlLoc, patchColPad);


bucknerIsoPatchCol(volidx==reconSubj,:) = []; 
bucknerDsPatchCol(volidx==reconSubj,:) = []; 
bucknerDsMaskPatchCol(volidx==reconSubj,:) = []; 

% load selected ADNI subject volumes
fnames = fullfile(SYNTHESIS_DATA_PATH, testdataset, 'md', [sys.usrname, '_restor_md_*']);
testmd = loadmd(fnames);
dsSubjNii = testmd.loadModality(dsSubjMod, reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(testmd.loadVolume(dsSubjMaskMod, reconSubj));
isoSubjInAtlNii = testmd.loadModality(isoSubjInAtlMod, reconSubj);
dsSubjInAtlNii = testmd.loadModality(dsSubjInAtlMod, reconSubj);
dsSubjInAtlMaskVol = testmd.loadVolume(dsSubjInAtlMaskMod, reconSubj);
subjInAtlTform = load(testmd.getModality(dsSubjInAtlMatMod, reconSubj));

% prepare necessary inputs for conditional-based reconstruction
subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
tform = subjInAtlTform.tform;
atlVolSize = size(dsSubjInAtlNii.img);
subjLoc2AtlSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims);
atlLoc2SubjSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims, 'backward');

%%

disp('knn quilt'); 

% crop the subvolumes from the iso, ds and mask volumes
isoSubjInAtlPatch = cropVolume(double(isoSubjInAtlNii.img), subvolLoc, subvolLoc + subvolSize - 1);
dsSubjInAtlPatch = cropVolume(double(dsSubjInAtlNii.img), subvolLoc, subvolLoc + subvolSize - 1);
dsSubjInAtlMaskPatch = cropVolume(double(dsSubjInAtlMaskVol), subvolLoc, subvolLoc + subvolSize - 1);


% divide the volumes into patch libraries
libPatches = patchlib.vol2lib(dsSubjInAtlPatch, atlPatchSize); 
libMasks = patchlib.vol2lib(dsSubjInAtlMaskPatch, atlPatchSize);


% find the matching patches 
dstfun = @(x, y) wtdst(x, y, atlPatchSize);
[pIdx, dist] = knnsearch([bucknerIsoPatchCol ones(size(bucknerIsoPatchCol))], [libPatches libMasks], 'Distance', dstfun); 
retrievedPatches = bucknerIsoPatchCol(pIdx, :); 

% quilt the patches
reconVol = patchlib.quilt(retrievedPatches, size(dsSubjInAtlPatch) - atlPatchSize + 1, atlPatchSize); 

% view the results
% linearlly interpolated subvolume, ground truth iso subvolume, reconstruced Volume
view3Dopt(dsSubjInAtlPatch, dsSubjInAtlMaskPatch, isoSubjInAtlPatch, reconVol)


%% compute isotropic gaussian mixture model

disp('gmm quilt'); 

% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);
gmmK = 15; 
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
keepr = 1; 

% compute the gaussian mixture model. TODO: should use wgmmfit with model0
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
gmmIso = fitgmdist(X, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'Options', gmmopt);
fprintf('gmm took %3.3f sec\n', toc);
gmmIso = wgmm.gmdist2wgmm(gmmIso);


extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);
dsSubjNii = testmd.loadModality(dsSubjMod, reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(testmd.loadVolume(dsSubjMaskMod, reconSubj));

[quiltedSubvolIso, reconLoc, cntvol] = papago.subvolRecon(gmmIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);




