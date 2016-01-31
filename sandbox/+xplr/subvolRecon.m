%% Test subvolume reconstruction
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
gmmK = 15; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 5; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'buckner';

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load *ground truth* data column from training
fnames = fullfile(SYNTHESIS_DATA_PATH, traindataset, 'md', [sys.usrname, '_restor_md_*']);
trainmd = loadmd(fnames);
[bucknerIsoPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(trainmd, 'brainIso2Ds5Us5sizeReg', atlPatchSize, atlLoc, patchColPad);

% load selected test subject volumes
fnames = fullfile(SYNTHESIS_DATA_PATH, testdataset, 'md', [sys.usrname, '_restor_md_*']);
testmd = loadmd(fnames);
dsSubjNii = testmd.loadModality('brainDs5Us5', reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(testmd.loadVolume('brainDs5Us5Mask', reconSubj));
dsSubjInAtlNii = testmd.loadModality('brainDs5Us5Reg', reconSubj);
dsSubjInAtlMaskVol = testmd.loadVolume('brainDs5Us5RegMask', reconSubj);
subjInAtlTform = load(testmd.getModality('brainDs5Us5RegMat', reconSubj));

% prepare necessary inputs for conditional-based reconstruction
subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
tform = subjInAtlTform.tform;
atlVolSize = size(dsSubjInAtlNii.img);
subjLoc2AtlSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims);
atlLoc2SubjSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims, 'backward');
extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);


%% compute gaussian mixture model
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
gmm = fitgmdist(X, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'Options', gmmopt);
fprintf('Gaussian mixture model took %3.3f sec\n', toc);
gmm = wgmm.gmdist2wgmm(gmm);

%% reconstruct patches in ADNI volume and quilt
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);
tic
[quiltedSubvol, reconLoc, cntvol] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);
toc

%% substitute and repeat
% TODO: need to substitute in atl. Or will that even make a difference? optimally, we're looking at
% atl-iso, and that hasn't been fantastic anyway.
quiltedSubvol
dsSubjVol
[quiltedSubvol, reconLoc, cntvol] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% test new reconstruction via "pre-computed" R.

R = cor2interpmat(size(atlLoc2SubjSpace{1}), subjLoc2AtlSpace); % move this to pre-computation (?)
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);
tic;
[quiltedSubvol_R, reconLoc_R, cntvol_R] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg, R);
toc;

%%

%% visualize
isoSubjVol = testmd.loadVolume('brainIso2Ds5Us5size', reconSubj);
cropIsoSubjVol = cropVolume(isoSubjVol, reconLoc, reconLoc + size(quiltedSubvol) - 1); 
cropIsoSubjVol(isnan(quiltedSubvol)) = nan;

dsSubjVolWNans = dsSubjVol;
cropDsSubjVolWNans = cropVolume(dsSubjVolWNans, reconLoc, reconLoc + size(quiltedSubvol) - 1); 
cropDsSubjVolWNans(isnan(quiltedSubvol)) = nan;

subjWeightVolWNans = dsSubjWeightVol*1;
cropSubjWeightVolWNans = cropVolume(subjWeightVolWNans, reconLoc, reconLoc + size(quiltedSubvol) - 1); 
cropSubjWeightVolWNans(isnan(quiltedSubvol)) = nan;

view3Dopt(cropIsoSubjVol, cropSubjWeightVolWNans, cropDsSubjVolWNans, quiltedSubvol, cntvol);

%% TODO: try higher keepr and mrf (use subject space patches but setup mrf on atlas space)