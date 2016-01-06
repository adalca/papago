%% setup
subspacesetup; 
patchSize = ones(1, 3) * 9;
location = LOC_LEFT_CORTEX; 
locpad = ones(1, 3) * 1;

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNIIs.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);
nandspatches = dspatches;
nandspatches(maskpatches < 0.75) = nan;

%% reconstruct nanpatches via PCA build from nanpatches via ALS algorithm
% compute PCA via ALS
[alleigvecs, allscores, allsigmas, t2, varexpl, mu] = pca(nandspatches, 'algorithm', 'als');

% extract top 95 eigenvectors
c = cumsum(varexpl./sum(varexpl));
f = find(c > 0.95, 1, 'first');

% reconstruct with PCA 
fakesigmas = allsigmas;
fakesigmas(f+1:end) = 0.0000000001;
reconvpatches = pcax.inpaint(nandspatches(5, :)', alleigvecs, fakesigmas)';

% visualize
figuresc(); spreadpatch = @(x) subspacetools.reshape3Dto2D(x, patchSize);
subplot(311); imagesc(spreadpatch(dspatches(5, :) + meanpatch), [0, 1]); title('correct');
subplot(312); imagesc(spreadpatch(nandspatches(5, :) + meanpatch), [0, 1]); title('deletions');
subplot(313); imagesc(spreadpatch(reconvpatches + meanpatch), [0, 1]); title('reonstructed');

%% reconstruct nanpatches via PCA build from isopatches via ALS algorithm
% compute PCA
[isoalleigvecs, isoallscores, isoallsigmas] = pca(isopatches);

% inpaint
patchno = 5;
reconvpatches = pcax.inpaint(nandspatches(patchno, :)', isoalleigvecs, isoallsigmas)';

% visualize
figuresc(); spreadpatch = @(x) subspacetools.reshape3Dto2D(x, patchSize);
subplot(311); imagesc(spreadpatch(dspatches(patchno, :)'), [-1, 1]); title('correct');
subplot(312); imagesc(spreadpatch(nandspatches(patchno, :)'), [-1, 1]); title('deletions');
subplot(313); imagesc(spreadpatch(reconvpatches), [-1, 1]); title('reonstructed');

%% impute via our iterative pcaimpute method
% parameters
K = 7;
nRep = 100;
[recovpatches, meanX] = pcaimpute(nandspatches, K, nRep, isopatches);

% visualize
patchno = 5;
figuresc(); spreadpatch = @(x) subspacetools.reshape3Dto2D(x, patchSize);
subplot(411); imagesc(spreadpatch(isopatches(patchno, :)'), [-1, 1]); title('correct');
subplot(412); imagesc(spreadpatch(nandspatches(patchno, :)'), [-1, 1]); title('deletions');
subplot(413); imagesc(spreadpatch(recovpatches(patchno, :)), [-1, 1]); title('reonstructed');
subplot(414); imagesc(spreadpatch(mean(isopatches)'), [-1, 1]); title('mean patch');
