% Experiment: gmm on iso data, get the cluster assignments, try to get
%     the same clusters from model3.

% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
% atlPatchSize = ones(1, 3) * 5; 
% atlLoc = [20, 25, 35];
gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 5; %1, 3
patchColPad = ones(1, 3) * 1;

% train and test datasets
traindataset = 'buckner';
testdataset = 'adni';

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


%% verify data
err = abs(bucknerIsoPatchCol - bucknerDsPatchCol);
plot(bucknerDsMaskPatchCol(:), err(:), '.');

l = linspace(0, 1, 100);
z = kernelRegress(bucknerDsMaskPatchCol(:), err(:), l, 0.05);
hold on; plot(l, z, '.-');

axislabels('patch weights', 'abs error', '');
legend({'data', 'kernel regression'});

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


%% get estimate from ds-wgmm

% get cluster assignment
probs = gmmIso.logpost(X, X*0+1);
[~, mi] = max(probs, [], 2);

% prepare data for w-gmm mstep
dsX = bsxfun(@minus, bucknerDsPatchCol, mean(bucknerDsPatchCol, 2));
dsW = bucknerDsMaskPatchCol + 0.00001;
gammank = zeros(size(X, 1), gmmK);
gammankidx = sub2ind(size(gammank), (1:size(X, 1))', mi(:));
gammank(gammankidx) = 1;


% TODO: recompute gmmIso with that cluster assignment.
gmmFakeIso = wgmm(gmmIso.mu*0, gmmIso.sigma*0, gmmIso.pi*0, gmmIso.sigmainv*0);
gmmFakeIso.covarUpdateMethod = 'model0';
gmmFakeIso.muUpdateMethod = 'model0';
gmmFakeIso.logpUpdateMethod = 'model0';
gmmFakeIso.sigmareg = 0.0000001;
gmmFakeIso.mstep(X, X*0+1, gmmK, gammank); % just do one update for the gmm. TODO: need witght

[quiltedSubvolFakeIso, reconLoc, cntvol] = papago.subvolRecon(gmmFakeIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% do w-gmm mstep
weightfact = 1e10;prod(patchColPad*2+1) * 15; % sigma-recon weight threshold
gmmIsoBasedDs = wgmm(gmmIso.mu*0, gmmIso.sigma*0, gmmIso.pi*0, gmmIso.sigmainv*0);
gmmIsoBasedDs.sigmaopt = weightfact;%, atlPatchSize};
gmmIsoBasedDs.sigmareg = 0.0000001;
% gmmIsoBasedDs.covarMergeMethod = 'freq-prior';
% gmmIsoBasedDs.covarMergeMethod = 'none';
gmmIsoBasedDs.mstep(dsX, dsW.^dsW, gmmK, gammank); % just do one update for the gmm. TODO: need witght

[quiltedSubvolIsoBasedDs, reconLoc, cntvol] = papago.subvolRecon(gmmIsoBasedDs, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% do w-gmm mstep
dsiX = bsxfun(@minus, bucknerDsInterpPatchCol, mean(bucknerDsInterpPatchCol, 2));
weightfact = 1e10; % sigma-recon weight threshold
gmmIsoBasedInterpDs = wgmm(gmmIso.mu*0, gmmIso.sigma*0, gmmIso.pi*0, gmmIso.sigmainv*0);
gmmIsoBasedInterpDs.sigmaopt = {weightfact, atlPatchSize};
gmmIsoBasedInterpDs.sigmareg = 0.0000001;
gmmIsoBasedInterpDs.covarMergeMethod = 'freq-prior';
% gmmIsoBasedDs.covarMergeMethod = 'none';
gmmIsoBasedInterpDs.mstep(dsiX, dsW, gmmK, gammank); % just do one update for the gmm. TODO: need witght

[quiltedSubvolIsoBasedInterpDs, reconLoc, cntvol] = papago.subvolRecon(gmmIsoBasedInterpDs, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% wtw computations
for i = 1:gmmK
    patchidx = (mi == i);
    subsW = dsW(patchidx, :);
    
%     subsW = (subsW > 0.95) + 0.000000000001;
    
    wtw(:,:,i) = subsW' * subsW;
    swtw(:,:,i) = wtw(:,:,i) ./ sqrt(size(patchidx, 1));
end

err = (abs(gmmFakeIso.sigma - gmmIsoBasedDs.sigma));
l = linspace(min(wtw(:)), max(wtw(:)), 100);
z = kernelRegress(wtw(:), err(:), l, 1);
scatter(wtw(:), err(:), '.'); hold on; plot(l, z);
axislabels('W''W', 'abs sigma err', 'sigma err at w''w locations');

% figure();
% hist(swtw(:), 100);

%% compute diffs

for i = 1:gmmK, 
    siso = gmmIso.sigma(:,:,i);
    sinterp = gmmIsoBasedInterpDs.sigma(:,:,i);
    
    df(i) = mean(msd(siso, sinterp));
    diff = abs(siso - sinterp);
    sc{i} = [siso, sinterp, diff]; % ./ (abs(siso + sinterp)/2)];
end

df = df ./ max(df);

view2D(sc, 'titles', arrayfunc(@(x) sprintf('%3.2f msd', x), df));
% view3Dopt(quiltedSubvolIso,quiltedSubvolFakeIso,quiltedSubvolIsoBasedDs, quiltedSubvolIsoBasedInterpDs);

%%
csel = [11, 1];
shortIso = wgmm(gmmIso.mu(csel, :), gmmIso.sigma(:,:, csel), gmmIso.pi(csel), gmmIso.sigmainv(:,:,csel));
shortFakeIso = wgmm(gmmIso.mu(csel, :), gmmFakeIso.sigma(:,:, csel), gmmFakeIso.pi(csel), gmmFakeIso.sigmainv(:,:,csel));
shortDsIso = wgmm(gmmIso.mu(csel, :), gmmIsoBasedDs.sigma(:,:, csel), gmmIsoBasedDs.pi(csel), gmmIsoBasedDs.sigmainv(:,:,csel));


[quiltedSubvolshortIso, reconLoc, cntvol] = papago.subvolRecon(shortIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvolshortFakeIso, reconLoc, cntvol] = papago.subvolRecon(shortFakeIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvolshortDsIso, reconLoc, cntvol] = papago.subvolRecon(shortDsIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% visualize results
% get real data subvolume
isoSubjVol = testmd.loadVolume('brainIso2Ds5Us2size', reconSubj);
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
    quiltedSubvolFakeIso, quiltedSubvolIsoBasedDs);


%%
z = [patchview.reshapeto2D(cat(2, gmmIso.sigma, zeros(size(X, 2), 100, gmmK))); 
patchview.reshapeto2D(cat(2, gmmIsoBasedDs.sigma, zeros(size(X, 2), 100, gmmK)))];
imagesc(z); colormap jet;
