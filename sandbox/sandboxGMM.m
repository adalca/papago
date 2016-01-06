%% setup
subspacesetup;
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_VENTRICLE_EDGE; 
locpad = o3 * 1;
maskthr = 0.75;

% error function
% TODO: compare different error functions
errfun = @(x, y) subspacetools.patcherror(x, y, @msd, 'subpatch', patchSize, o3*3, o3*7);
errfun = @(x, y) subspacetools.patcherror(x, y, @msd, 'gauss', patchSize, 1.5);

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNIIs.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);
nandspatches = dspatches;
nandspatches(maskpatches < maskthr) = nan;
nNans = sum(maskpatches(:) < maskthr);
fprintf('%d / %d are nan (%3.2f%%)\n', nNans, numel(maskpatches), nNans./numel(maskpatches));

%% Gaussian Mixture and Visualization
nClust = 5;
X = isopatches; %(251:end, :);
patchElem = prod(patchSize);
patchSizeIn2D = [patchSize(1), patchSize(2)*patchSize(3)];

% fit a GMM.
% st = statset('Display', 'iter', 'MaxIter', 50);
% gmd = gmdistribution.fit(X, nClust, 'Replicates', 5, 'Regularize', 0.00001, 'Options', st);
[~, ci] = max(gmd.posterior(X), [], 2); % get cluster assignment

% sample from Gaussian mixtures and display
samples = gmmSampleGrid(gmd, patchSizeIn2D, [10, 2]);

% PCA of individual clusters
[coeff, scores, latent, expl, mu, sampscores, oppscores] = deal(cell(nClust, 1));
for c = 1:nClust
    cidx = ci == c;
    [coeff{c}, scores{c}, latent{c}, ~, expl{c}, mu{c}] = pca(X(cidx, :));
    
    % check usage of pcax.project.
    cdata = bsxfun(@minus, X(cidx, :)', mu{c}');
    scorescheck = pcax.project(cdata, coeff{c})';
    assert(max(max(abs(scorescheck - scores{c}))) < 1e-6);
    
    % project all opposite points % TODO: why using different mu?
    cnotdata = bsxfun(@minus, X(~cidx, :)', mu{c}');
    oppscores{c} = pcax.project(cnotdata, coeff{c});
end

% sample from Gaussian mixtures using only first two components, and display
pcafun = @(coeff, latent) coeff * (diag([latent(1:2); zeros(numel(latent) - 2, 1)])) * coeff';
pcaCov = cellfunc(pcafun, coeff, latent); 
[samplesCropCov, samplesCropCovCell] = gmmSampleGrid(gmd, patchSizeIn2D, [10, 2], cat(3, pcaCov{:}));

% project all points in this cluster
for c = 1:nClust
    sampscores{c} = pcax.project(bsxfun(@minus, samplesCropCovCell{c}, gmd.mu(c,:))', coeff{c})';
end

%% Visualization
% displaying clusters, statistics, samples
nRows = 5;
nClustShow = 5;
[compProp, cps] = sort(gmd.ComponentProportion, 'descend');
showclusts = cps(1:nClustShow);

% display means and covariances
figuresc(); colormap gray; 
for i = 1:numel(showclusts)
    c = showclusts(i);
    subplot(nRows, nClustShow, i); 
    muim = reshape(gmd.mu(c, :), patchSizeIn2D);
    imagesc(muim); axis equal; axis off; 
    title(sprintf('cluster %d (%f%%) mean', c, compProp(i)));
    
    subplot(nRows, nClustShow, nClustShow + i); 
    imagesc(reshape(gmd.Sigma(:,:,c), [patchElem, patchElem])); axis equal; axis off; 
    title(sprintf('cluster %d covar', c));
end

% show samples from clusters
for i = 1:numel(showclusts)
    c = showclusts(i);
    subplot(nRows, nClustShow, nClustShow*2 + i); 
    imagesc(samples{c}(:,:,1)); axis equal;  axis off; colorbar; 
    title('Sample Patches using full Covariances'); 
end

% PCA of individual clusters
for i = 1:numel(showclusts)
    c = showclusts(i);
    cidx = ci == c;
    
    subplot(nRows, nClustShow, nClustShow*3 + i); 
    plot(scores{c}(:, 1), scores{c}(:, 2), '.'); hold on;
    plot(oppscores{c}(:, 1), oppscores{c}(:, 2), 'o'); axis equal;
    title(sprintf('PCA space for cluster %d (2 comp expl %3.2f%%)', c, sum(expl{c}(1:2))));
    legend({sprintf('cluster %d (%d)', c, sum(cidx)), sprintf('not cluster %d (%d)', c, sum(~cidx))});
end

% display sampling from PCA
for i = 1:numel(showclusts)
    c = showclusts(i);
    
    % project points onto PCA and draw them on PCA plots
    subplot(nRows, nClustShow, nClustShow*3 + i); 
    plot(sampscores{c}(:,1), sampscores{c}(:,2), '*'); 
    legend({sprintf('cluster %d (%d)', c, sum(cidx)), ...
        sprintf('not cluster %d (%d)', c, sum(~cidx)), ...
        sprintf('Samples from PCA space')});
        
    % show samples as matrix
    subplot(nRows, nClustShow, nClustShow*4 + i); 
    imagesc(samplesCropCov{c}(:,:,1)); axis equal;  axis off; colorbar;
    title('Sample Patches using First 2 Princ. Cmp'); 
end

%% try overall pca for display purposes.
[ocoeff, oscores, ~, ~, oexpl, omu] = pca(X);
figuresc(); plot(oscores(:, 1), oscores(:, 2), '.');

%% try to project (without pca) from hidden data.
Xi = gmmInpaint(nandspatches, gmd);
gmmResultExamples(nandspatches, dspatches, Xi, isopatches, gmd, errfun, patchSize, ci, 2);

%% try gmm from isotropic and then gmm impute with isotropic nan'ed.
nClust = 50;
nanisopatches = isopatches; nanisopatches(maskpatches < maskthr) = nan;
st = statset('Display', 'iter', 'MaxIter', 50);
gmd = gmdistribution.fit(isopatches, nClust, 'Replicates', 5, 'Regularize', 0.00001, 'Options', st);
[~, ci] = max(gmd.posterior(isopatches), [], 2); % get cluster assignment

% impute without pca.
Xi = gmmInpaint(nanisopatches, gmd);
gmmResultExamples(nandspatches, dspatches, Xi, isopatches, gmd, errfun, patchSize, ci, 2);

% TODO: so, with good initialization we're likely to do well. So maybe do higher scale and work way
% down? That way can have bigger patches? and the mrf/voting.
blurisovol = cellfunc(@(x) volblur(x.img, 2), niis.iso);
blurisopatches = subspacetools.vols2patchcol(blurisovol,  patchSize, location, locpad);
bnanisopatches = blurisopatches; bnanisopatches(~maskpatches) = nan;

% impute without pca.
Xi = gmmInpaint(bnanisopatches, gmd);
gmmResultExamples(nandspatches, dspatches, Xi, isopatches, gmd, errfun, patchSize, ci, 2);
% should also see blurisopatches

% TODO ideas: 
% -- do some sort of learn in low res but with high res info...? n/s. Warning: this is not
%   enforcing the structure you know!
% -- is there a way to learn with *no* info and project with the blurry info? N/S
% -- test via edges (if edges are well detected)
%
% -- instead of project in probability, do exemplar project to a series of exemplars. but we don't
% have exemplars. but parts of them somehow?

%% test: blur iso patches, and in ~missing put high res. 
%   - Train gmm and impute (1 step) (equiv to gmmimpute with nRep = 1?
%   - do gmmimpute

% TODO/MUST DO: instead of using hard maps, use probability maps and use those when you do cluster
% assignments?

% train data
blurisovol = cellfunc(@(x) volblur(x.img, 2), niis.iso);
blurisopatches = subspacetools.vols2patchcol(blurisovol,  patchSize, location, locpad);
bnanisopatches = blurisopatches; bnanisopatches(maskpatches < maskthr) = nan;
bhrisopatches = blurisopatches; bhrisopatches(maskpatches > maskthr) = isopatches(maskpatches > maskthr);

opts = {'valinit', bhrisopatches, 'trueData', isopatches, 'realErrFun', errfun};

% Test 1: learn new gmm and impute
[Xi, gmd] = gmmimpute(nanisopatches, nClust, 0, 'maxIter', 1, opts{:});
fprintf('%e --> %e\n', errfun(bhrisopatches, isopatches), errfun(Xi, isopatches))

% compute some distances
% note: cheating with ci?
gmmResultExamples(nanisopatches, bhrisopatches, Xi, isopatches, gmd, errfun, patchSize, ci, 2);

% Test 2: full gmm+impute several iterations
[Xi, gmd] = gmmimpute(nanisopatches, nClust, 5, 'maxIter', 100, opts{:});
fprintf('%e --> %e\n', errfun(bhrisopatches, isopatches), errfun(Xi, isopatches));

[Xiw, gmdw] = gmmimpute(nanisopatches, nClust, 5, 'maxIter', 100, opts{:}, 'valwts', maskpatches);
fprintf('%e --> %e\n', errfun(bhrisopatches, isopatches), errfun(Xiw, isopatches));
% TODO: including weights and using wt method does worse, but:
% 1. it may be better when we need to actually use dspatches, and don't cheat (using isopatches)! [Done: YES! though not better by too much]
% 2. maybe don't use wt for posteriors, just for assignment? [Done: NOPE]

% conclusion: method seems to work! so hope to go hierarchical

% TODO: think of katie's suggestion of gmm with rotation.

