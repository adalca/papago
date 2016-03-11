
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
atlLoc = [52,43,34];


atlPatchSize = ones(1, 3) * 6; 
atlLoc = [69, 51, 43]; % scale 3
atlLoc = [53, 48, 48];

atlPatchSize = ones(1, 3) * 5;
atlLoc = [20, 25, 35];
% atlLoc = [22, 20, 20]*2;

gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 7; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'adni';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

subvolSize = atlPatchSize + 5;

% modalities
ds = 5;
us = 2;
isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
dsSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Reg', ds, us);
dsSubjInAtlMaskMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2MaskReg', ds, us);
dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);

dsSubjMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2', ds, us);
dsSubjMaskMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Mask', ds, us);
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);


% isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
% dsSubjInAtlMod = sprintf('brainDs%dUs%dReg', ds, us);
% dsSubjInAtlMaskMod = sprintf('brainDs%dUs%dRegMask', ds, us);
% dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
% dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);
% 
% dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
% dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);
% isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);

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
[~, mi] = max(isopost, [], 2);


%% 
threshold = 0.5;
clear sigmast meanst sigmas0 means0 sigmas means means_iso sigmas_iso
for k = 1:gmmK
    dfn = @(w) diag((w<0.5)*1000000000);
    dfn = @(w) diag(-log(w+eps)*100000000);
    dfn = @(w) diag(0.0001*(exp(-10*(w-1)) - 0.99));

    % get isotropic mu, sigma
    X = bucknerIsoPatchCol(mi == k, :);
    X = bsxfun(@minus, X, mean(X, 2));
    sigmast(:,:,k) = cov(X);
    meanst(k, :) = mean(X);
        
    % get downsampled mu, sigma
    X0 = bucknerDsPatchCol(mi == k, :);
    X0 = bsxfun(@minus, X0, mean(X0, 2));
    sigmas0(:,:,k) = cov(X0);
    means0(k, :) = mean(X0);
    
    % get weights
    W = bucknerDsMaskPatchCol(mi == k, :);
    
    tic;
    wgs{k} = wgmmfit(X0, W, 1, 'model4fn', dfn,  'replicates', 1, 'initmethod', 'model4exp05', 'regularizationValue', 0, 'TolFun', -inf, 'MaxIter', 5);
    fprintf('wgmm cluster %d done in %5.3fs\n', k, toc);
%     
%     tic;
%     wgs_iso{k} = wgmmfit(X, W, 1, 'model4fn', dfn,  'replicates', 1, 'initmethod', 'model4exp05', 'regularizationValue', 0, 'TolFun', -inf, 'MaxIter', 5);
%     fprintf('wgmm cluster %d done in %5.3fs\n', k, toc);
    
    % set up downsampled data with nans in any region with W < 0.75
    X0(W<threshold) = nan;
    % estimate mean and sigmas with ECM algorithm.
    tic;
    [means(k, :), sigmas(:,:,k)] = ecmnmlex(X0, 'twostage', 5, 0.000); 
    fprintf('cluster %d done in %5.3fs\n', k, toc);
    
%     tic;
%     X(W<threshold) = nan;
%     [means_iso(k, :), sigmas_iso(:,:,k)] = ecmnmlex(X, 'twostage', 5, 0.000); 
%     fprintf('cluster %d done in %5.3fs\n', k, toc);
    
%     imagesc([sigmas_iso(:,:,k), wgs{k}.sigma(:,:,1), sigmas_iso(:,:,k) - wgs{k}.sigma(:,:,1)])
%     drawnow;

%     % set up downsampled data with nans in any region with W < 0.75
%     X0tmp = X0; 
%     Xtmp = X; 
%     X0tmp(W<threshold) = nan;
%     Xtmp(W<threshold) = nan; 
%     
%     
%     hackNum = 50; 
%     [Wvals, idxWvals] = sort(W, 'descend'); 
%     nanColumns = sum(~isnan(X0tmp),1) < hackNum; 
%     colVals = find(nanColumns) ; 
%     inds = sub2ind( size(X0tmp),  idxWvals(1:hackNum, nanColumns), repmat(colVals, [hackNum, 1]) ); 
%     X0tmp( inds(:) ) = X0( inds(:) ); 
%     Xtmp( inds(:) ) = X( inds(:) ); 
%     
%     [means(k, :), sigmas(:,:,k)] = ecmnmlex(X0tmp, 'twostage', 10, 0.001);
end

for k = 1:gmmK
    wgmm_m4ds_means(k, :) = wgs{k}.mu;
    wgmm_m4ds_sigmas(:, :, k) = wgs{k}.sigma;
    
    wgmm_m4iso_means(k, :) = wgs{k}.mu;
    wgmm_m4iso_sigmas(:, :, k) = wgs{k}.sigma;
end

%% setup wgmms
wgt = wgmm(meanst, sigmast + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
wg0 = wgmm(means0, sigmas0 + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
wg = wgmm(means, sigmas + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
wg_iso = wgmm(means_iso, sigmas_iso + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
wg_m4ds = wgmm(wgmm_m4ds_means, wgmm_m4ds_sigmas + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
wg_m4iso = wgmm(wgmm_m4iso_means, wgmm_m4iso_sigmas + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));

%% 

% reconstruct 
[quiltedSubvolt_k5, reconLoc, cntvol] = papago.subvolRecon(wgt, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvol0_k5, reconLoc, cntvol] = papago.subvolRecon(wg0, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvol_k5, reconLoc, cntvol] = papago.subvolRecon(wg, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvol_iso_k5, reconLoc, cntvol] = papago.subvolRecon(wg_iso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvol_m4ds_k5, reconLoc, cntvol] = papago.subvolRecon(wg_m4ds, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

[quiltedSubvol_m4iso_k5, reconLoc, cntvol] = papago.subvolRecon(wg_m4iso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

view3Dopt(quiltedSubvolt_k5, quiltedSubvol0_k5, quiltedSubvol_iso_k5, quiltedSubvol_k5, quiltedSubvol_m4ds_k5, quiltedSubvol_m4iso_k5)

%%
qsel = quiltedSubvolt_k5;
isoSubjVol = testmd.loadVolume(isoSubjMod, reconSubj);
cropIsoSubjVol = cropVolume(isoSubjVol, reconLoc, reconLoc + size(qsel) - 1); 
cropIsoSubjVol(isnan(qsel)) = nan;

% get linearly-interpolated data
dsSubjVolWNans = dsSubjVol;
cropDsSubjVolWNans = cropVolume(dsSubjVolWNans, reconLoc, reconLoc + size(qsel) - 1); 
cropDsSubjVolWNans(isnan(qsel)) = nan;

% get subject weight
subjWeightVolWNans = dsSubjWeightVol*1;
cropSubjWeightVolWNans = cropVolume(subjWeightVolWNans, reconLoc, reconLoc + size(qsel) - 1); 
cropSubjWeightVolWNans(isnan(qsel)) = nan;

view3Dopt(cropIsoSubjVol, cropDsSubjVolWNans, cropSubjWeightVolWNans, quiltedSubvolt_k5, quiltedSubvol0_k5, quiltedSubvol_k5)


