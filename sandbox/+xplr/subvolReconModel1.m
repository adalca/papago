% Experiment: gmm on iso data, get the cluster assignments, try to get
%     the same clusters from model3.

% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
atlPatchSize = ones(1, 3) * 7; 
atlLoc = [20, 25, 35];
gmmK = 15; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 5; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'buckner';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

% subvolSize = atlPatchSize + 5;

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
gmmIsoBasedDs.sigmaopt = {weightfact, atlPatchSize};
gmmIsoBasedDs.sigmareg = 0.0000001;
gmmIsoBasedDs.covarMergeMethod = 'freq-prior';
% gmmIsoBasedDs.covarMergeMethod = 'none';
gmmIsoBasedDs.mstep(dsX, dsW, gmmK, gammank); % just do one update for the gmm. TODO: need witght

[quiltedSubvolIsoBasedDs, reconLoc, cntvol] = papago.subvolRecon(gmmIsoBasedDs, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


%% do w-gmm mstep with model1
weightfact = 1e10;prod(patchColPad*2+1) * 15; % sigma-recon weight threshold
gmmIsoBasedDsM1 = wgmm(gmmIso.mu*0, gmmIso.sigma*0, gmmIso.pi*0, gmmIso.sigmainv*0);
gmmIsoBasedDsM1.sigmaopt = weightfact;%, atlPatchSize};
gmmIsoBasedDsM1.sigmareg = 0.0000001;
gmmIsoBasedDsM1.covarUpdateMethod = 'memsafe-model1';
% gmmIsoBasedDsM1.covarMergeMethod = 'freq-prior';
% gmmIsoBasedDsM1.covarMergeMethod = 'none';
gmmIsoBasedDsM1.mstep(dsX, dsW, gmmK, gammank); % just do one update for the gmm. TODO: need witght

clear wtw;
for i = 1:gmmK
    patchidx = (mi == i);
    subsW = dsW(patchidx, :);
    
%     subsW = (subsW > 0.95) + 0.000000000001;
    
    wtw(:,:,i) = subsW' * subsW;
end





[quilted, ~, ~, modReconLocs, reconPatches] = papago.subvolRecon(gmmIsoBasedDsM1, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);
% 
% gmmIsoBasedDsW1 = wgmm(gmmIso.mu*0, wtw, gmmIsoBasedDsM1.pi);
% warning('need to select logp properly!');
% [pweights, reconLoc, cntvol, a, b] = papago.subvolRecon(gmmIsoBasedDsW1, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
%     dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjWeightVol*1, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);
% 


%%

% subvolSize = atlPatchSize;
clear err bw

for sid = 1:10
    reconSubj = sid;
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

    [quilted, reconLoc, ~, modReconLocs, reconPatches] = papago.subvolRecon(gmmIsoBasedDsM1, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
        dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

    [quiltedSubvolIso, reconLoc, cntvol] = papago.subvolRecon(gmmIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
        dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

    
    % get real data subvolume
    isoSubjVol = testmd.loadVolume('brainIso2Ds5Us2size', reconSubj);
    cropIsoSubjVol = cropVolume(isoSubjVol, reconLoc, reconLoc + size(quilted) - 1); 
    cropIsoSubjVol(isnan(quilted)) = nan;

    % get subject weight
    subjWeightVolWNans = dsSubjWeightVol*1;
    cropSubjWeightVolWNans = cropVolume(subjWeightVolWNans, reconLoc, reconLoc + size(quilted) - 1); 
    cropSubjWeightVolWNans(isnan(quilted)) = nan;
    cropSubjWeightVol4bw = cropSubjWeightVolWNans;
    cropSubjWeightVol4bw(isnan(cropSubjWeightVol4bw)) = 0;

    truep{sid} = cropIsoSubjVol;
    reconp{sid} = quilted;
    qiso{sid} = quiltedSubvolIso;

    err{sid} = abs(quilted - quiltedSubvolIso);
    bw{sid} = bwdist(cropSubjWeightVol4bw);
    bwn{sid} = bwdist(isnan(cropSubjWeightVolWNans));
    bwn{sid}(cropSubjWeightVolWNans==1) = 0;
    bw{sid}(isnan(cropSubjWeightVolWNans)) = nan;
    sid
end

errv = cellfunc(@(x) x(:), err);
bwv = cellfunc(@(x) x(:), bw);
bwnv = cellfunc(@(x) x(:), bwn);
allerr = cat(1, errv{:});
allbw= cat(1, bwv{:});
allbwn= cat(1, bwnv{:});
