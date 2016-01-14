%% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;

gmmmethod = 'groundtruth'; % 'compute', 'load', 'groundtruth'
method = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward

subject = 3; %1, 3

%% load volumes
% we'll use some local, small subset, of the buckner data for now.

% TODO: should move this to using the existing md structure, 
% once both katie and adrian have access to the same large data.
bpref = sprintf('buckner%02d', subject);

% subject space
subjisonii = loadNii(sprintf('%s/%s/%s/%s_brain_cropped5.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % iso
dsnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds 
dsmasknii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_dsmask.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));

% atlas space (i.e. registered data)
isoregnii = loadNii(sprintf('%s/%s/%s/%s_brain_iso_2_ds5_us5_size_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds
dsregnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds
tf = load(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_reg.mat', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));
dsmaskregnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_dsmask_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));
    
%% prepare necessary inputs for reconstruction
% atlLoc, atlPatchSize, subjVol, subjWeightVol, atlLoc2SubjSpace, method
subjVol = double(dsnii.img);
subjWeightVol = logical(dsmasknii.img);

subjDims = dsnii.hdr.dime.pixdim(2:4);
atlDims = dsregnii.hdr.dime.pixdim(2:4);
subjLoc2AtlSpace = tform2cor3d(tf.tform, size(subjVol), subjDims, size(dsregnii.img), atlDims);
atlLoc2SubjSpace = tform2cor3d(tf.tform, size(subjVol), subjDims, size(dsregnii.img), atlDims, 'backward');
extraReconArg = ifelse(strcmp(method, 'inverse'), subjLoc2AtlSpace, regVal);

%% "groundtruth" testing
% atlMu, atlSigma: done below

% get the covariance from the (isotropic) patch we're looking at in atlas space and see if we can
% reconstruction the right patch in subject space.
isoSubjVolInAtl = double(isoregnii.img);

% get mu and sigma in atlas space
atlPatch = cropVolume(isoSubjVolInAtl, atlLoc, atlLoc + atlPatchSize - 1);
atlMu = zeros(atlPatchSize);
atlSigma = atlPatch(:) * atlPatch(:)' + eye(prod(atlPatchSize)) * regVal;

% reconstruct
reconPatch_gt = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
    subjVol, subjWeightVol, atlLoc2SubjSpace, method, extraReconArg);

%% load data column
% load *ground truth* data column
atlasesDataSetup; % atlases
BUCKNER_PATH_ORIG = fullfile(BUCKNER_PATH, 'orig'); % original files
BUCKNER_PATH_PROC = fullfile(SYNTHESIS_DATA_PATH, 'buckner/proc'); % buckner processing 
md = restorationmd(dsAmounts, BUCKNER_PATH_PROC, SYNTHESIS_DATA_PATH, 'buckner');
[col, layeridx, volidx] = subspacetools.md2patchcol(md, 'brainIso2Ds5Us5size', atlPatchSize, atlLoc, [2, 2, 2]);
col(volidx == subject, :) = [];

%% single - gaussian estimation
gmm = fitgmdist(col, 1);
atlMu = gmm.mu(:)';
atlSigma = gmm.Sigma + eye(prod(atlPatchSize)) * regVal;

reconPatch_sg = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
    subjVol, subjWeightVol, atlLoc2SubjSpace, method, extraReconArg);

%% gaussian mixture model
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
% minus this subject.
K = 5;
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model
tic
gmm = fitgmdist(col, K, 'regularizationValue', regVal, 'replicates', 10, 'Options', gmmopt);
fprintf('Gaussian mixture model took %3.3f sec\n', toc);

% reconstruct
% TODO: really, we should only be computing logp, not the reconstruction!
logps = zeros(1, K);
reconPatches = cell(1, K);
subjPatchLocs = cell(1, K);
for k = 1:K
    atlMu = gmm.mu(k, :)';
    atlSigma = gmm.Sigma(:, :, k) + eye(prod(atlPatchSize)) * regVal;
    
    [reconPatches{k}, logps(k), subjPatchLocs{k}] = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
        subjVol, subjWeightVol, atlLoc2SubjSpace, method, extraReconArg);
end

% get the optimal patch reconstruction
[~, clust] = max(log(gmm.ComponentProportion) + logps);
reconPatch_kg = reconPatches{clust};
subjPatchLoc = subjPatchLocs{clust};

%% visualize and compare
subjPatchSize = size(reconPatch_gt);

% compare subject patches for: true isotropic volume, linear interpolation, mask.
dsSubjPatch = cropVolume(dsnii.img, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
dsSubjPatch(isnan(reconPatch_gt)) = nan;
maskSubjPatch = cropVolume(dsmasknii.img, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
maskSubjPatch(isnan(reconPatch_gt)) = nan;
correctSubjPatch = cropVolume(subjisonii.img, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
correctSubjPatch(isnan(reconPatch_gt)) = nan;

% 3D visualization
view3Dopt(correctSubjPatch, maskSubjPatch, dsSubjPatch, ...
    reconPatch_gt, reconPatch_sg, reconPatch_kg);

% 2D visualization
figuresc();
subplot(2, 3, 1); imagesc(patchview.reshapeto2D(correctSubjPatch)); title('true');
subplot(2, 3, 2); imagesc(patchview.reshapeto2D(maskSubjPatch)); title('mask');
subplot(2, 3, 3); imagesc(patchview.reshapeto2D(dsSubjPatch));  title('linear interp');
subplot(2, 3, 4); imagesc(patchview.reshapeto2D(reconPatch_gt));   title('our recon groundtruth');
subplot(2, 3, 5); imagesc(patchview.reshapeto2D(reconPatch_sg));   title('our recon single gaussian');
subplot(2, 3, 6); imagesc(patchview.reshapeto2D(reconPatch_kg));   title('our recon gmm');
colormap gray;

%% TODO: compare logp in this manner with others.
