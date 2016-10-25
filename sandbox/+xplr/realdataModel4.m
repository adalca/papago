% Real data test of model4-wgmm

% initialize
setup

%% parameters
nSimSamples = 4750;

atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
%atlPatchSize = ones(1, 3) * 5; 
%atlLoc = [20, 25, 35];
gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 7; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'buckner';
testdataset = 'buckner';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

subvolSize = atlPatchSize + 5;

% modalities
ds = 5;
us = 5;
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
[bucknerIsoPatchCol, layeridx, volidx] = ...
    subspacetools.md2patchcol(trainmd, isoSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsMaskPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMaskMod, atlPatchSize, atlLoc, patchColPad);


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

%% learn a gmm
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model. TODO: should use wgmmfit with model0
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
gmdistIso = fitgmdist(X, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'Options', gmmopt);
fprintf('gmm took %3.3f sec\n', toc);
gmmIso = wgmm.gmdist2wgmm(gmdistIso);

isopost = gmdistIso.posterior(X);

%% real ds data.
[n, d] = size(bucknerDsMaskPatchCol);
r = randsample(n, nSimSamples);
Xt = bsxfun(@minus, bucknerIsoPatchCol(r, :), mean(bucknerIsoPatchCol(r, :), 2));
X0 = bsxfun(@minus, bucknerDsPatchCol(r, :), mean(bucknerDsPatchCol(r, :), 2));
W = bucknerDsMaskPatchCol(r, :);
W = max(W, 0.000001);

fn3 = @(w) diag((-log(w))) * Adj * diag((-log(w))) * 0.0006;

% wg = wgmm.fit(X0, W, gmmK, 'replicates', 2, 'model4fn', fn3); % ! with defaults set to model 4.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3, 'updateMethod', 'model5', 'model4fn', fn3);
[~, gammank, Xhat1] = wg.estep(X0, W);
xhaterr1 = msd(Xt, Xhat1, 2);


% see how the estep behaves given the right mstep
wg2 = wgmm(gmmIso.mu, gmmIso.sigma, gmmIso.pi); % ! with defaults set to model 4.
wg2.model4fn = fn3;
[~, ~, Xhat2] = wg2.estep(X0, W);
xhaterr2 = msd(Xt, Xhat2, 2);

wg3 = wgmm(gmmIso.mu*0, gmmIso.sigma*0, gmmIso.pi*0); % ! with defaults set to model 4.
wg3.sigmareg = 0.000001;
wg3.mstep(Xt, W, gmmK, isopost(r, :)); % re-estimate mu, sigma, pi given right Xt and assignment (this should be perfect).
wg3.model4fn = fn3;
[~, ~, Xhat3] = wg3.estep(X0, W);
xhaterr3 = msd(Xt, Xhat3, 2);

figure(); 
l = linspace(0, 0.015, 50); 
c0 = hist(msd(Xt, X0, 2), l); plot(l, c0); hold on; title('histogram of errors');
c1 = hist(xhaterr1, l); plot(l, c1); 
c2 = hist(xhaterr2, l); plot(l, c2); 
c3 = hist(xhaterr3, l); plot(l, c3); 
legend({sprintf('linear %3.4f,', mean(msd(Xt, X0, 2))), ...
    sprintf('model4 %3.4f', mean(xhaterr1)), ...
    sprintf('m4-estep from gmmIso params %3.4f', mean(xhaterr2)), ...
    sprintf('m4-mstep from Xt and isoposteriors %3.4f', mean(xhaterr3))});

%% reconstruct subvolume
[quiltedSubvolM4, reconLoc, cntvol] = papago.subvolRecon(wg, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

%% visualize results

% get real data subvolume
isoSubjVol = testmd.loadVolume('brainIso2Ds5Us2size', reconSubj);
cropIsoSubjVol = cropVolume(isoSubjVol, reconLoc, reconLoc + size(quiltedSubvolM4) - 1); 
cropIsoSubjVol(isnan(quiltedSubvolM4)) = nan;

% get linearly-interpolated data
dsSubjVolWNans = dsSubjVol;
cropDsSubjVolWNans = cropVolume(dsSubjVolWNans, reconLoc, reconLoc + size(quiltedSubvolM4) - 1); 
cropDsSubjVolWNans(isnan(quiltedSubvolM4)) = nan;

% get subject weight
subjWeightVolWNans = dsSubjWeightVol*1;
cropSubjWeightVolWNans = cropVolume(subjWeightVolWNans, reconLoc, reconLoc + size(quiltedSubvolM4) - 1); 
cropSubjWeightVolWNans(isnan(quiltedSubvolM4)) = nan;

view3Dopt(cropIsoSubjVol, cropSubjWeightVolWNans, cropDsSubjVolWNans, quiltedSubvolM4)

%%
i = randsample(nSimSamples, 1);
imagesc(patchview.reshapeto2D([Xt(i, :); X0(i, :); Xhat(i, :);  Xhat2(i, :)], atlPatchSize), [0, 0.5]); colormap gray;

view3Dopt(reshape(Xt(i, :), atlPatchSize), reshape(X0(i, :), atlPatchSize), ...
    reshape(W(i, :), atlPatchSize), reshape(Xhat(i, :), atlPatchSize));


%% compute model3 wgmm from linearly-interpolated  data with weights
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
gmmDs = wgmmfit(X, W, gmmK, 'regularizationValue', regVal, 'replicates', 3, wgmmopts{:}, 'updateMethod', 'model3');
fprintf('Ds wgmm took %3.3f sec\n', toc);

[quiltedSubvolDs] = papago.subvolRecon(gmmDs, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);
