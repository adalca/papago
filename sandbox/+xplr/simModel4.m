% Simulate data and test aspects of model4-wgmm

% initialize
setup

%% parameters
nSimSamples = 1000;

atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
atlPatchSize = ones(1, 3) * 5; 
atlLoc = [20, 25, 35];
gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 7; %1, 3
patchColPad = ones(1, 3) * 1;

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

%% simulate destorying some of the data according to the model

[Xt, cl] = gmmIso.sample(nSimSamples);
W = rand(size(Xt));
W(W > 0.3) = 1;
W(W < 0.1) = 0.1;

X0 = zeros(size(Xt));
for i = 1:size(Xt, 1)
    w = W(i, :);
    D = diag((-log(w)).^2);
    X0(i, :) = mvnrnd(Xt(i, :), D);
end

x0err = msd(Xt(:), X0(:));

%% E-step (mixture of wiener filters) check

% try just making a new wgmm, and do estep with the correct sigma, mu
wg = wgmm(gmmIso.mu, gmmIso.sigma, gmmIso.pi); % ! with defaults set to model 4.
[~, ~, Xhat] = wg.estep(X0, W);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));

%% run full model4-wgmm on this simulated data.
% estimate gmm.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3); % ! with defaults set to model 4.

% estimate Xhat and cluster assignments
[~, gammank, Xhat] = wg.estep(X0, W);
[~, cidxhat] = max(gammank, [], 2);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));

figure();
hist(cidxhat, 1:gmmK);

% TODO: use wgmm model 3 to get sigma, mu, and weiner filter to get xhat, and compare

%% simulate data from cluster, but with real weights
[Xt, cl] = gmmIso.sample(nSimSamples);

W = bucknerDsMaskPatchCol(randsample(size(bucknerDsMaskPatchCol, 1), nSimSamples), :);
W = max(W, 0.000001);

X0 = zeros(size(Xt));
for i = 1:size(Xt, 1)
    w = W(i, :);
    D = diag((-log(w)).^2);
    X0(i, :) = mvnrnd(Xt(i, :), D);
end

% estimate gmm.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3); % ! with defaults set to model 4.

% estimate Xhat and cluster assignments
[~, gammank, Xhat] = wg.estep(X0, W);
[~, cidxhat] = max(gammank, [], 2);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));


%% new D function using adjacency
q = size2ndgridvec(atlPatchSize);
Adj = exp(-pdist2(q, q));
fn1 = @(w) diag((-log(w)).^2);
fn2 = @(w) Adj * diag((-log(w)).^2);
fn3 = @(w) diag((-log(w))) * Adj * diag((-log(w))) * 0.0006;
fn4 = @(w) (-log(w))' * (-log(w)) .* Adj * 0.0006;
% view2D({D1, D2; D3 D4}, 'titles', {'D=diag(-log(W)).^2', 'Adj*D', 'sqrt(D)*Adj*sqrt(D)', '(W''W)*Adj'}); colormap gray;

W = bucknerDsMaskPatchCol(randsample(size(bucknerDsMaskPatchCol, 1), nSimSamples), :);
W = max(W, 0.000001);

X0 = zeros(size(Xt));
for i = 1:size(Xt, 1)
    w = W(i, :);
    %D = diag((-log(w)).^2);
    D = fn3(w);
    X0(i, :) = mvnrnd(Xt(i, :), D);
end

% estimate gmm.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3, 'model4fn', fn3); % ! with defaults set to model 4.

% estimate Xhat and cluster assignments
[~, gammank, Xhat] = wg.estep(X0, W);
[~, cidxhat] = max(gammank, [], 2);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));


