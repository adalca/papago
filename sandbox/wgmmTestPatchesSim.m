%% setup
subspacesetup; 
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_VENTRICLE_EDGE; LOC_LEFT_CORTEX;
locpad = o3 * 2;

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);

%% simulation using patches
% this way we'll have some patch-like-looking data
N = 1000;                   % number of patches to sample from GMM
D = size(isopatches, 2);    % dimensionality of the patch
K = 3;                      % number of clusters in GMM to learn from data

% learn and sample from a gmm from the nLearn random patches.
nlearn = min(size(isopatches, 1), max(10000, D)+1);
r = randperm(size(isopatches, 1));
Xlearn = isopatches(r(1:nlearn), :);
[X, mu, sigma, pi, cidx, sigmainv] = subspacetools.simgmm(N, D, K, 'frompatchesgmm', Xlearn);

% put the sampled patches along with weights into weighted gmm class
W = ones(size(X)); 
wgmmsampled = wgmm(X, W, K, mu, sigma, pi, sigmainv); 

% visualize the means of the learned GMM 
figure(); imagesc(subspacetools.reshapeN3Dto2D(wgmmsampled.mu, patchSize)); colormap gray; title('means'); axis equal off;

% visualize some new samples vs real patches
subspacetools.comparePatchGMMData(Xlearn, ones(nlearn, 1), X, ones(N, 1), patchSize, 30);
subplot(1,2,1); title('real data'); subplot(1,2,2); title('resampled data');

%% run wgmm with weights = 1 on simulation using patches
W = ones(size(X));
figure();
[isosimfwgmm] = wgmmfit(X, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);

% visualize resulting isofwgmm and compare with gmm from which we sampled.
isosimfwgmm.visualize(wgmmsampled);
figure(); imagesc(subspacetools.reshapeN3Dto2D(isosimfwgmm.mu, patchSize)); colormap gray; title('means'); axis equal off;
[sX, sidx] = isosimfwgmm.sample(N);
subspacetools.comparePatchGMMData(X, cidx(:), sX, sidx(:), patchSize, 10);
subplot(K,2,1); title('input data'); subplot(K,2,2); title('wgmm resample data');

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

%% try reconstructions
prtile = 95;
reconeigPatches = papago.recon(isofwgmm, isopatches, W, 'eig', prtile);
% TODO: figure out how to do the reconstruction so that you only invert coeffs once, not every time
% invert coeffs * W
reconweigPatches = papago.recon(isofwgmm, isopatches, W, 'weig', prtile); 
reconpcaPatches = papago.recon(isofwgmm, isopatches, W, 'pca', prtile);

% visualize
rX = {reconeigPatches, reconeigPatches, reconeigPatches};
papago.visRecon(isopatches, W, rX, patchSize, 'titles', {'EIG', 'WEIG', 'PCA'});

%%
W = ones(size(X));
[isosimfwgmm] = wgmmfit(X, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);
W = max(maskpatches(1:N, :), 0.001);
[isosimfwgmm2] = wgmmfit(X, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);
W = max(maskpatches(1:N, :), 0.001);
[isosimfwgmm3] = wgmmfit(X./W, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);
[isosimfwgmm4] = wgmmfit(X.*W, W, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 5);

gmmfit = fitgmdist(X, K, 'RegularizationValue', 0.00001, 'replicates', 5);
gmfit.sigma = gmmfit.Sigma;

%%
gmms = {wgmmsampled, isosimfwgmm, isosimfwgmm2, isosimfwgmm3, isosimfwgmm4, isofwgmm9};
titles = {'sample sigma', 'W=1 wgmm', 'usual Ws wgmm', 'usual Ws, wgmm with X/W', 'usual Ws, wgmm with X*W', 'qq'}; 

gmms = {wgmmsampled, isosimfwgmm, isosimfwgmm2, isofwgmm9, gmfit};
titles = {'sample sigma', 'W=1 wgmm', 'usual Ws wgmm', 'qq', 'W=1 via fitgmdist'}; 

T = numel(gmms);
for i = 1:K; 
    kstr = sprintf(' sigma for mixture %d', i);
    ktitles = cellfunc(@(s) [s, kstr], titles);
    for t = 1:T
        subplot(T, K, sub2ind([K, T], i, t)); 
        imagesc(gmms{t}.sigma(:,:,i)); title(ktitles{t});
    end
end
colormap gray;
