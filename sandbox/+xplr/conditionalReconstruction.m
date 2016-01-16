%% Test Conditional Reconstruction methods
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
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
subjVol = double(dsSubjNii.img);
subjWeightVol = logical(adnimd.loadVolume('brainDs5Us5Mask', reconSubj));
dsSubjInAtlNii = adnimd.loadModality('brainDs5Us5Reg', reconSubj);
subjInAtlTform = load(adnimd.getModality('brainDs5Us5RegMat', reconSubj));

% prepare necessary inputs for conditional-based reconstruction
subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
tform = subjInAtlTform.tform;
atlVolSize = size(dsSubjInAtlNii.img);
subjLoc2AtlSpace = tform2cor3d(tform, size(subjVol), subjDims, atlVolSize, atlDims);
atlLoc2SubjSpace = tform2cor3d(tform, size(subjVol), subjDims, atlVolSize, atlDims, 'backward');
extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);

%% "groundtruth" testing
% Test the ability to reconstruct a patch, where we compute a (regularized) covariance from the
% single "true" (isotropic) subject patch that's hidden from the reconstruction method.

% get "true" mu and sigma in atlas space
isoSubjInAtlVol = adnimd.loadVolume('brainIso2Ds5Us5sizeReg', reconSubj);
atlPatch = cropVolume(isoSubjInAtlVol, atlLoc, atlLoc + atlPatchSize - 1);
atlMu = zeros(atlPatchSize);
atlSigma = atlPatch(:) * atlPatch(:)' + eye(prod(atlPatchSize)) * regVal;

% reconstruct the given patch
reconPatch_gt = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
    subjVol, subjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArg);

%% single - gaussian estimation
gmm = fitgmdist(bucknerIsoPatchCol, 1);
atlMu = gmm.mu(:)';
atlSigma = gmm.Sigma + eye(prod(atlPatchSize)) * regVal;

reconPatch_sg = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, ...
    subjVol, subjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArg);

%% gaussian mixture model
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model
tic
gmm = fitgmdist(bucknerIsoPatchCol, gmmK, 'regularizationValue', regVal, 'replicates', 10, 'Options', gmmopt);
fprintf('Gaussian mixture model took %3.3f sec\n', toc);

% reconstruct
% TODO: really, we should only be computing logp, not the reconstruction!
logps = zeros(1, gmmK);
logp2 = zeros(1, gmmK);
reconPatches = cell(1, gmmK);
subjPatchLocs = cell(1, gmmK);
for k = 1:gmmK
    atlMu = gmm.mu(k, :)';
    atlSigma = gmm.Sigma(:, :, k) + eye(prod(atlPatchSize)) * regVal;
    
    [reconPatches{k}, subjPatchLocs{k}, logps(k)] = paffine.recon(atlMu, atlSigma, atlLoc, ...
        atlPatchSize, subjVol, subjWeightVol, atlLoc2SubjSpace, crmethod, extraReconArg);
end

% get the optimal patch reconstruction
[~, clust] = max(log(gmm.ComponentProportion) + logps);
reconPatch_kg = reconPatches{clust};

%% get some useful subject-space data for comparison or initialization
% get the subject space coordinates
[~, subjPatchLoc, subjPatchSize] = ...
    paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);

% compare subject patches for: true isotropic volume, linear interpolation, mask.
dsSubjPatch = cropVolume(subjVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
dsSubjPatch(isnan(reconPatch_gt)) = nan;
maskSubjPatch = cropVolume(subjWeightVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
maskSubjPatch(isnan(reconPatch_gt)) = false;

isoSubjVol = adnimd.loadVolume('brainCropped5', reconSubj);
correctSubjPatch = cropVolume(isoSubjVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
correctSubjPatch(isnan(reconPatch_gt)) = nan;

%% visualize and compare

% 3D visualization
view3Dopt(correctSubjPatch, maskSubjPatch*1, dsSubjPatch, ...
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

%% 3D visualize of the gmm options, ssd to "true" data, and logp based on sparse data.
fprintf('logP  :');
fprintf('%9.2f ', logps);
fprintf('\nlogPv2:');
fprintf('%9.2f ', logp2);
fprintf('\nssd   :');
valid = ~isnan(correctSubjPatch);
ssds = cellfun(@(x) ssd(x(valid), correctSubjPatch(valid)), reconPatches);
fprintf('%9.2f ', ssds);
newline;

view3Dopt(correctSubjPatch, reconPatches{:})
