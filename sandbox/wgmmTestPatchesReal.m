%% setup
subspacesetup; 
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_LEFT_CORTEX; LOC_VENTRICLE_EDGE
locpad = o3 * 2;

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);

%% WGMM on real patches isotropic data with weights = 1
K = 7;
W = ones(size(isopatches));
figure();
[isofwgmm] = wgmmfit(isopatches, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);
isofwgmm.visualize();
figure(); imagesc(subspacetools.reshapeN3Dto2D(isofwgmm.mu, patchSize)); colormap gray; title('means'); axis equal off;

% visualize new samples from resulting isowgmm.
[sX, sidx] = isofwgmm.sample(N);
subspacetools.comparePatchGMMData(isopatches, randi([1, K], size(isopatches, 1), 1), sX, sidx, patchSize, 5);
subplot(K,2,1); title('input data'); subplot(K,2,2); title('wgmm resample data');

%% WGMM on real patches *isotropic* data with real weights
K = 7;
W = min(1, maskpatches + 0.001);
figure();
[isofwgmm] = wgmmfit(isopatches, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);
isofwgmm.visualize();
figure(); imagesc(subspacetools.reshapeN3Dto2D(isofwgmm.mu, patchSize)); colormap gray; title('means'); axis equal off;

% visualize new samples from resulting isowgmm.
[sX, sidx] = isofwgmm.sample(N);
subspacetools.comparePatchGMMData(isopatches, randi([1, K], size(isopatches, 1), 1), sX, sidx, patchSize, 5);
subplot(K,2,1); title('input data'); subplot(K,2,2); title('wgmm resample data');

%% run WGMM on ds patchs
K = 3;

% blurvols = cellfunc(@(x) volblur(x.img, 1.5), niis.iso);
[bisopatches, bisopatchidx] = subspacetools.vols2patchcol(blurvols, patchSize, location, locpad);
dspatches = maskpatches .* isopatches + (1 - maskpatches) .* bisopatches;

% TODO must do! Try X = X_o .* (1/W). and see how covariance compares with normal gmm of X_o.
W = min(1, maskpatches + 0.001);
figure();
[dsfwgmm] = wgmmfit(dspatches, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 3);
dsfwgmm.visualize();
figure(); imagesc(subspacetools.reshapeN3Dto2D(dsfwgmm.mu, patchSize)); colormap gray; title('means'); axis equal off;

% visualize new samples from resulting isowgmm.
[sX, sidx] = dsfwgmm.sample(N);
subspacetools.comparePatchGMMData(dspatches, randi([1, K], size(dspatches, 1), 1), sX, sidx, patchSize, 5);
subplot(K,2,1); title('input data'); subplot(K,2,2); title('wgmm resample data');

%% reconstuct via pca projection
% here doing eig of sigma, rather than pca of data. Maybe this is partly why this is better?
% TODO: try some sort of weighted PCA projection?
% TODO: weight this patch controbution to whole volume based on loglik.
prtile = 95;
reconeigPatches = papago.recon(dsfwgmm, dspatches, W, 'eig', prtile);
% reconweigPatches = papago.recon(dsfwgmm, dspatches, W, 'weig', prtile);
reconpcaPatches = papago.recon(dsfwgmm, dspatches, W, 'pca', prtile);

% visualize
rX = {dspatches, reconeigPatches, reconeigPatches, reconeigPatches};
papago.visRecon(isopatches, W, rX, patchSize, 'titles', {'DS', 'EIG', 'WEIG', 'PCA'});

%% try: pca beforehand.
% [coeff, score, latent, q, r, mu] = pca(isopatches);
% X = isopatches;
% W = X * 0 + 1;
% [isofwgmm] = wgmmfit(X, W, K, 'debug', false, 'replicates', 5);

%% need to use NN instead of linear. 
% look into bucknerDataProcess do NN via registeredVolumeInterp in register();
