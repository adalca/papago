%% compare different models in reconstructing subvolumes.
% TODOs: 
%   - the experiment of calling gmm on iso, then get the cluster assignments, then try to get
%     the same clusters from model3.

% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlPatchSize = ones(1, 3) * 5; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
atlLoc = [20, 23, 36];
gmmK = 15; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 3; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'adni';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

% modalities
ds = 5;
us = 2;
isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
dsSubjInAtlMod = sprintf('brainDs%dUs%dReg', ds, us);
dsSubjInAtlMaskMod = sprintf('brainDs%dUs%dRegMask', ds, us);
dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);

dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);


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
% [bucknerDsInterpPatchCol, ~, ~] = ... 
%     subspacetools.md2patchcol(trainmd, dsInterpSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);

% load selected ADNI subject volumes
fnames = fullfile(SYNTHESIS_DATA_PATH, testdataset, 'md', [sys.usrname, '_restor_md_*']);
testmd = loadmd(fnames);
dsSubjNii = testmd.loadModality(dsSubjMod, reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(testmd.loadVolume(dsSubjMaskMod, reconSubj));
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
extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);

%% compute isotropic gaussian mixture model
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model. TODO: should use wgmmfit with model0
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
gmmIso = fitgmdist(X, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'Options', gmmopt);
fprintf('gmm took %3.3f sec\n', toc);
gmmIso = wgmm.gmdist2wgmm(gmmIso);

[quiltedSubvolIso, reconLoc, cntvol] = papago.subvolRecon(gmmIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% compute wgmm from linearly-interpolated  data with weights
weightfact = prod(patchColPad*2+1) * 15; % sigma-recon weight threshold
weightfact = size(bucknerDsPatchCol, 1) ./ 5;
% weightfact = 1000;
wgmmopts = {'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', {weightfact, atlPatchSize}, 'verbose', 1, 'covarMergeMethod', 'freq-prior'};
%wgmmopts = {'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', {weightfact}, 'verbose', 1, 'covarMergeMethod', 'wfact-mult-adapt'};

% compute the gaussian mixture model
warning('todo: add spatial smoothness prior to sigma? Read Eugenio Paper!');
tic
X = bsxfun(@minus, bucknerDsPatchCol, mean(bucknerDsPatchCol, 2));
W = bucknerDsMaskPatchCol + 0.0000000001;
W = W.^3;
% W(W < 0.9) = 0.0000000000001;
gmmDs = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, wgmmopts{:});
fprintf('Ds wgmm took %3.3f sec\n', toc);

[quiltedSubvolDs] = papago.subvolRecon(gmmDs, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% compute wgmm from true (iso) data with weights
weightfact = prod(patchColPad*2+1) * 15; % sigma-recon weight threshold

% compute the gaussian mixture model
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
W = bucknerDsMaskPatchCol + 0.0000000001;
W(W < 0.7) = 0.0000000000001;
gmmIsoW = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', weightfact, 'verbose', 1);
fprintf('Iso W wgmm took %3.3f sec\n', toc);

[quiltedSubvolIsoW, reconLoc, cntvol] = papago.subvolRecon(gmmIsoW, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% compute wgmm from a mix of true (iso) data (in high-weight areas) and linearly-interp data.
tic
W = bucknerDsMaskPatchCol;
X = bucknerIsoPatchCol;
X(W < 0.1) = bucknerDsPatchCol(W < 0.1); % ds data in low-areas.
X = bsxfun(@minus, X, mean(X, 2));
W = W + 0.00001;
gmmIsoDsW = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', weightfact, 'verbose', 1);
fprintf('Iso W gaussian mixture model took %3.3f sec\n', toc);

[quiltedSubvolIsoDsW, reconLoc, cntvol] = papago.subvolRecon(gmmIsoDsW, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% compute wgmm from a mix of true (iso) data (in high-weight areas) and interp-based rotated data.

warning('we need to look at these patches. Currently, ds-interp has 0s in between interpolated areas');
tic
W = bucknerDsMaskPatchCol;
X = bucknerDsInterpPatchCol;
X(W < 0.1) = bucknerDsPatchCol(W < 0.1); % ds data in low-areas.
X = bsxfun(@minus, X, mean(X, 2));
W = W + 0.00001;
gmmDsinterpDsW = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', weightfact, 'verbose', 1);
fprintf('Iso W gaussian mixture model took %3.3f sec\n', toc);

[quiltedSubvolDsinterpDsW, reconLoc, cntvol] = papago.subvolRecon(gmmDsinterpDsW, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% compute wgmm from linearly-interpolated  data with high-weights
tic
X = bsxfun(@minus, bucknerDsPatchCol, mean(bucknerDsPatchCol, 2));
W = bucknerDsMaskPatchCol;
W(W < 0.9) = 0.00001;
gmmDsThrW = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'MaxIter', 20, 'TolFun', 0.001, 'regularizationWeight', weightfact, 'verbose', 1);
fprintf('Iso W gaussian mixture model took %3.3f sec\n', toc);

[quiltedSubvolDsThrW] = papago.subvolRecon(gmmDsThrW, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% visualize results
% get real data subvolume
isoSubjVol = testmd.loadVolume(isoSubjMod, reconSubj);
cropIsoSubjVol = cropVolume(isoSubjVol, reconLoc, reconLoc + size(quiltedSubvolIso) - 1); 
cropIsoSubjVol(isnan(quiltedSubvolIso)) = nan;

% get linearly-interpolated data
dsSubjVolWNans = dsSubjVol;
cropDsSubjVolWNans = cropVolume(dsSubjVolWNans, reconLoc, reconLoc + size(quiltedSubvolIso) - 1); 
cropDsSubjVolWNans(isnan(quiltedSubvolIso)) = nan;

% get subject weight
subjWeightVolWNans = dsSubjWeightVol*1;
cropSubjWeightVolWNans = cropVolume(subjWeightVolWNans, reconLoc, reconLoc + size(quiltedSubvolIso) - 1); 
cropSubjWeightVolWNans(isnan(quiltedSubvolIso)) = nan;


view3Dopt(cropIsoSubjVol, cropSubjWeightVolWNans, cropDsSubjVolWNans, quiltedSubvolIso, ...
    quiltedSubvolDs, quiltedSubvolIsoW, quiltedSubvolIsoDsW, quiltedSubvolDsinterpDsW, ...
    quiltedSubvolDsThrW, cntvol);
