%% Test subvolume reconstruction
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE-10; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 3; %1, 3
patchColPad = ones(1, 3) * 2;

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load *ground truth* data column from buckner
bucknermd = loadmd([SYNTHESIS_DATA_PATH, 'buckner', '_restor_md_*']);
[bucknerIsoPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(bucknermd, 'brainIso2Ds5Us5size', atlPatchSize, atlLoc, patchColPad);

% load selected ADNI subject volumes
adnimd = loadmd([SYNTHESIS_DATA_PATH, 'adni', '_restor_md_*']);
dsSubjNii = adnimd.loadModality('brainDs5Us5', reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(adnimd.loadVolume('brainDs5Us5Mask', reconSubj));
dsSubjInAtlNii = adnimd.loadModality('brainDs5Us5Reg', reconSubj);
dsSubjInAtlMaskVol = adnimd.loadVolume('brainDs5Us5RegMask', reconSubj);
subjInAtlTform = load(adnimd.getModality('brainDs5Us5RegMat', reconSubj));

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
gmm = fitgmdist(bucknerIsoPatchCol, gmmK, 'regularizationValue', regVal, 'replicates', 10, 'Options', gmmopt);
fprintf('Gaussian mixture model took %3.3f sec\n', toc);

%% reconstruct patches in ADNI volume and quilt
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);
[quiltedSubvol, reconLoc, cntvol] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% visualize
isoSubjVol = adnimd.loadVolume('brainIso2Ds5Us5size', reconSubj);
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