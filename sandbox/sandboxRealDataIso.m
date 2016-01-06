%% setup
subspacesetup; 
patchSize = ones(1, 3) * 5;
location = [56, 72, 92];
locpad = ones(1, 3) * 1;
nDestroy = 75;

%% prepare data (patches)
if ~exist('niis', 'var'); load('data/bucknerNIIs.mat'); end
[patches, patchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);

%% run PCA on high res patches
[alleigvecs, allscores, allsigmas, ~, varexpl] = pca(patches(1:4000, :));

% extract top 95 eigenvectors
c = cumsum(varexpl./sum(varexpl));
f = find(c > 0.99, 1, 'first');

eigenvecs = alleigvecs(:, 1:f);
sigmas = allsigmas(1:f);
scores = allscores(:, 1:f);

%% destroy and reconstruct
nReps = 3;

% Maybe just substitute last sigmas?
warning('should not reconstruct with full eigenvecs and eigenvals');
sandboxDeleteReconstruct(patches, scores, alleigvecs, allsigmas, patchSize, nDestroy, nReps);

%% reconstruct the PCA space
% destroy data in an inconsistent way
nanpatches = subspacetools.destroyFeatures(patches, nDestroy, 'rand-inconsistent');
testpatches = nanpatches(end-99:end, :);
[alleigvecs, allscores, allsigmas, t2, varexpl] = pca(nanpatches, 'algorithm', 'als');

% extract top 95 eigenvectors
c = cumsum(varexpl./sum(varexpl));
f = find(c > 0.95, 1, 'first');

eigenvecs = alleigvecs(:, 1:f);
sigmas = allsigmas(1:f);
scores = allscores(1:f);

%% destroy and reconstruct with new PCA for last 100 samples
nDestroy = 75;
nReps = 1;

fakesigmas = allsigmas;
fakesigmas(f+1:end) = 0;
fakesigmas = fakesigmas + 0.0000000001;
[recovpatches, nanpatches] = sandboxDeleteReconstruct(patches, allscores, alleigvecs, fakesigmas, patchSize, nDestroy, nReps);

% visualize
figuresc()
patchno = 5;
figuresc(); spreadpatch = @(x) subspacetools.reshape3Dto2D(x, patchSize);
subplot(215); imagesc(spreadpatch(patches(patchno, :)'), [-1, 1]); title('correct');
subplot(216); imagesc(spreadpatch(nanpatches(patchno, :)'), [-1, 1]); title('deletions');
subplot(217); imagesc(spreadpatch(recovpatches(patchno, :)), [-1, 1]); title('reonstructed');
subplot(218); imagesc(spreadpatch(mean(isopatches)'), [-1, 1]); title('mean patch');

%% reconstruct with pcaimpute.
K = 7;
nRep = 100;
nanpatches = subspacetools.destroyFeatures(patches, nDestroy, 'rand-inconsistent');
[recovpatches{1}, meanX] = pcaimpute(nanpatches, K, nRep, patches);

% visualize 
figuresc()
patchno = 5;
figuresc(); spreadpatch = @(x) subspacetools.reshape3Dto2D(x, patchSize);
subplot(215); imagesc(spreadpatch(patches(patchno, :)'), [-1, 1]); title('correct');
subplot(216); imagesc(spreadpatch(nanpatches(patchno, :)'), [-1, 1]); title('deletions');
subplot(217); imagesc(spreadpatch(recovpatches(patchno, :)), [-1, 1]); title('reonstructed');
subplot(218); imagesc(spreadpatch(mean(isopatches)'), [-1, 1]); title('mean patch');
