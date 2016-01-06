%% setup
subspacesetup; 
nComp = 5;                  % number of model components
patchSize = [5, 5];         % patch size (dimentionality)
assert(nComp <= prod(patchSize)); 
nSamples = 1000;            % number of simulated samples
nDestroy = round(prod(patchSize) * 4/5);
nReps = 3;
invalidvox = [true(1, nDestroy), false(1, prod(patchSize) - nDestroy)];

%% simulate the model (PCA eigenvalues and eigenvectors)

% get an orthonormal basis.
alleigvecs = dct(eye(prod(patchSize)))'; 
alleigvecs = [alleigvecs(:, 2:nComp + 1), alleigvecs(:, [1, nComp + 2:end])]; % move DC to end

% get the actual non-zero basis & eigenvalues.
allsigmas = [rand(1, nComp), zeros(1, prod(patchSize) - nComp)];
sigmas = allsigmas(1:nComp);
eigenvecs = alleigvecs(:, 1:nComp); 

% compute the intrinsic space covariance: U * diag(Sigma) * U^(-1)
covar = alleigvecs * diag(allsigmas) / alleigvecs; 
% verify that the covariance makes sense
[verU, verSigma] = eig(covar); 
assert(all(isclose(sort(allsigmas, 'descend')', real(diag(verSigma)))));

%% Create sample data
scores = mvnrnd(zeros(1, nComp), sigmas * eye(nComp), nSamples);
patches = scores * eigenvecs';

% Sanity check: recover the score for this patch in one of several ways
recoveredScores1 = pcax.project(patches', eigenvecs)';
assert(all(all(isclose(recoveredScores1, scores))));
% (P(25x1)*E^T)*inv(E(25x5)*E^T) = S(1x5) |||||    X * U * (U'*U)^(-1)
recoveredScores2 = (patches * eigenvecs) / (eigenvecs'*eigenvecs); 
assert(all(all(isclose(recoveredScores2, scores))));

%% destroy and recover data
sandboxDeleteReconstruct(patches, alleigvecs, allsigmas, patchSize, nDestroy, nReps);

%% recover pca space from *consistently* deleted features

% destroy in a consistent way
nanpatches = subspacetools.destroyFeatures(patches, nDestroy, 'rand-consistent');

% build PCA space via ALS
[coeff, recovscore, latent] = pca(nanpatches, 'algorithm', 'als');

% sort eigenvalues
[srtsigmas, si] = sort(allsigmas, 'descend');
[srtrecovsigmas, sri] = sort(latent, 'descend');

% visualize results
figure();
subplot(231); imagesc(alleigvecs(:, si)); title('original U');
subplot(232); imagesc(real(coeff(:, sri))); title('estimated U');
subplot(233); plot(real(srtrecovsigmas), '-.'); hold on; plot(srtsigmas, '-*'); title('Eigenvals');
subplot(234); imagesc(scores(:, si(1:5))); title('original scores');
subplot(235); imagesc(real(recovscore(:, sri(1:5)))); title('recovered scores');

%% recover pca space from *inconsistently* deleted features, and recover features

% destroy data in an inconsistent way
nanpatches = subspacetools.destroyFeatures(patches, nDestroy, 'rand-inconsistent');

% build PCA space via ALS
[coeff, recovscore, latent] = pca(nanpatches, 'algorithm','als');
[srtsigmas, si] = sort(allsigmas, 'descend');
[srtrecovsigmas, sri] = sort(latent, 'descend');

% impute via our iterative pcaimpute method
K = 5;
[reconpcaimpute, meanpatch] = pcaimpute(nanpatches, K, 100, patches);

% recover with simple pca reconstruction
reconpca = pcax.recon(coeff(:, 1:end), recovscore(:, 1:end)', meanpatch)';

% visualize results
figure();
subplot(331); imagesc(alleigvecs(:, si)); title('original U');
subplot(332); imagesc(real(coeff(:, sri))); title('estimated U');
subplot(333); plot(real(srtrecovsigmas), '-.'); hold on; plot(srtsigmas, '-*'); title('Eigenvals');
subplot(334); imagesc(scores(:, si(1:5))); title('original scores');
subplot(335); imagesc(real(recovscore(:, sri(1:5)))); title('recovered scores');
subplot(336); imagesc(reshape(patches(patchno, :)', patchSize), [-1, 1]); title('correct');
subplot(337); imagesc(reshape(nanpatches(patchno, :)', patchSize), [-1, 1]); title('deletions');
subplot(338); imagesc(reshape(reconpcaimpute(patchno, :), patchSize), [-1, 1]); title('pcaimpute-recon');
subplot(338); imagesc(reshape(reconpcaimpute(patchno, :), patchSize), [-1, 1]); title('pca-recon');
