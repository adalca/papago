ed%% Explore doing a PCA decomposition and GMM as an initialization to a full GMM.
%% Setup
o3 = ones(1, 3);
patchSize = o3 * 7;
locpad = o3 * 1;
location = ceil(rand(1, 3) .* (size(mfptrs{1}, 'rigidRegBrain') - patchSize - 10));
K = 25;

% gmm parameters
gmmopt = statset('Display', 'final', 'MaxIter', 10);
gmmargs = {'RegularizationValue', regval, 'replicates', 5, 'Options', gmmopt};

minExplained = 95;
% TODO: also try synthetic mm data (?)



%% load data
tic;
mfptrs = arrayfunc(@(x) matfile(md.getModality('matfile', x)), 1:md.getNumSubjects);
colcmd = @(v, loc, pad) subspacetools.matfile2patchcol(mfptrs, v, patchSize, loc, pad);
% [isopatches, layeridx, volidx] = colcmd('rigidRegBrain', location, locpad);
fprintf('Loading data took %3.2fs\n', toc);

%% normal gmm
tic;
gmdistfull = fitgmdist(isopatches, K, gmmargs{:}); 
[~, gmonlyll] = gmdistfull.posterior(isopatches);
gmonlytoc = toc;

%% pca, gmm on pca, then full gmm
% initialize full gmm based on posterior assignments
tic; 
[~, score, ~, ~, explained] = pca(isopatches);
cexplained = cumsum(explained);
f = find(cexplained > minExplained, 1, 'first');
fprintf('%d components explain %3.1f%% of variance\n', f, cexplained(f));

% compute gmm dist on smaller dimensions
pcagmdist = fitgmdist(score(:, 1:f), K, gmmargs{:}); 

% use this gmm to initialize a full gmm. Note, we can only run one replicate.
% TODO: perhaps one idea is to run several replicates by doing slight perturbations of the posterior
post = pcagmdist.posterior(score(:, 1:f));
[~, mi] = max(post, [], 2);

% for those clusters that have zero votes, give them some reasonable points
h = hist(mi, 1:K); 
z = find(h == 0);
allbigclusters = find(h > (numel(z)+1));
ptsinbigclusters = find(any(bsxfun(@eq, mi, allbigclusters), 2));
for i = 1:numel(z), 
    k = z(i); 
    [~, pti] = max(post(ptsinbigclusters, k)); 
    mi(ptsinbigclusters(pti)) = k; 
    ptsinbigclusters(pti) = [];
end
h = hist(mi, 1:K); assert(all(h > 0));

gmdistAfterPCA1 = fitgmdist(isopatches, K, gmmargs{:}, 'Start', mi, 'replicates', 1);
[~, pca1ll] = gmdistAfterPCA1.posterior(isopatches);
pca1toc = toc;


%% pca, gmm on pca, then full gmm
% initialize full gmm based on mu and sigma estimates
tic; 
[coeff, score, latent, tquared, explained] = pca(isopatches);
cexplained = cumsum(explained);
f = find(cexplained > minExplained, 1, 'first');
fprintf('%d components explain %3.1f%% of variance\n', f, cexplained(f));
gmdist = fitgmdist(score(:, 1:f), K, gmmargs{:}); 
q = gmdist.posterior(score(:, 1:f));

% use this gmdist to initlize mu and sigma.
mu = zeros(K, size(isopatches, 2));
sigma = zeros(size(isopatches, 2), size(isopatches, 2), K);
for k = 1:K
    mu(k, :) = q(:, k)' * isopatches ./ sum(q(:, k));
    d = bsxfun(@minus, isopatches, mu(k, :));
    dp = bsxfun(@times, q(:, k), d);
    sigma(:,:,k) = dp' * d ./ sum(q(:, k)) + eye(size(isopatches, 2)) * 0.0001;
end

% compute full gmm
S = struct('mu', mu, 'Sigma', sigma, 'ComponentProportion', gmdist.ComponentProportion);
gmdistAfterPCA2 = fitgmdist(isopatches, K, gmmargs{:}, 'Start', S, 'replicates', 1);
[~, pca2ll] = gmdistAfterPCA2.posterior(isopatches);
pca2toc = toc;

%% analyse
fprintf('%15s%15s%15s%15s\n', 'Stat', 'directGMM', 'PCA-init-post', 'PCA-init-norm');
fprintf('%15s%14.3fs%14.3fs%14.3fs\n', 'Time Taken', gmonlytoc, pca1toc, pca2toc);
fprintf('%15s%15.3f%15.3f%15.3f\n', 'NLogLik', gmonlyll, pca1ll, pca2ll);

% visualize normals and sigmas

c2 = oogreedymatching(gmdistfull.mu, gmdistAfterPCA1.mu);
c3 = oogreedymatching(gmdistfull.mu, gmdistAfterPCA2.mu);

figuresc(); colormap gray; 
dispK = min(K, 5);
subplot(3, dispK+1, 1); imagesc(gmdistfull.mu); 
title('direct-GMM mus');
subplot(3, dispK+1, dispK+2); imagesc(gmdistAfterPCA1.mu(c2, :)); 
title('pca-GMM init via pt assignment. mus');
subplot(3, dispK+1, 2*(dispK+1) + 1); imagesc(gmdistAfterPCA2.mu(c3, :));  
title('pca-GMM init via mu,sigma. mus');
for k = 1:dispK
    ttl = sprintf('direct-GMM K:%d', k);
    subplot(3, dispK+1, 1+k); imagesc(gmdistfull.Sigma(:,:,k)); 
    axis equal; axis off; title(ttl);
    ttl = sprintf('pca-GMM init via pt assignment. K:%d', k);
    subplot(3, dispK+1, dispK+2+k); imagesc(gmdistAfterPCA1.Sigma(:,:,c2(k))); 
    axis equal; axis off; title(ttl);
    ttl = sprintf('pca-GMM init via mu,sigma. K:%d', k);
    subplot(3, dispK+1, 2*(dispK+1)+1+k); imagesc(gmdistAfterPCA2.Sigma(:,:,c3(k))); 
    axis equal; axis off; title(ttl);
end

% ss: conclusion: 
% second method is the fastest, significantly faster in some cases depending on K, repetitions
% results are not very stable, **but** the second and third method can even give higher loglik than
% the original run. Conclusion? try the pca-init! saves time and gives good results
