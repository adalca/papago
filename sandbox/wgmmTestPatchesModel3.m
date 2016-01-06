%% setup
% wgmmTestPatchesModel3.m
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_VENTRICLE_EDGE; 
locpad = o3 * 2;

% extract patches
if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);

blurvols = cellfunc(@(x) volblur(x.img, 2), niis.iso);
[bisopatches, bisopatchidx] = subspacetools.vols2patchcol(blurvols, patchSize, location, locpad);

%% test model 3 compared to matlab's normal gmmfit, data from patches
K = 3;
X = isopatches;
isogmmfit = fitgmdist(X, K, 'RegularizationValue', 0.00001, 'replicates', 10);
isogmfit = wgmm(X, maskpatches, K, isogmmfit.mu, isogmmfit.Sigma, isogmmfit.ComponentProportion);


% wgmm with W with extra division
W = max(maskpatches, 0.01);
test1args = {'covarUpdateMethod', 'model3', 'muUpdateMethod', 'model3', 'logpUpdateMethod', 'model3'};
[isowgmmm3] = wgmmfit(X, W, K, 'debug', false, test1args{:}, 'sigmareg',  0.0001, 'replicates', 20);


dspatches = W .* isopatches + (1 - W) .* bisopatches;
X = dspatches;
dsgmmfit = fitgmdist(X, K, 'RegularizationValue', 0.00001, 'replicates', 10);
dsgmfit = wgmm(X, maskpatches, K, dsgmmfit.mu, dsgmmfit.Sigma, dsgmmfit.ComponentProportion);

% wgmm with W with extra division
test1args = {'covarUpdateMethod', 'model3', 'muUpdateMethod', 'model3', 'logpUpdateMethod', 'model3'};
[dswgmmm3] = wgmmfit(X, W, K, 'debug', false, test1args{:}, 'sigmareg',  0.0001, 'replicates', 10);

%%
gmms = {isogmfit, isowgmmm3, dsgmfit, dswgmmm3};
titles = {'iso fitgmdist', 'iso realW model3', 'ds fitgmdist', 'ds m3'}; 

% show sigmas and mus
figuresc(1); figuresc(2); 
T = numel(gmms);
for i = 1:K; 
    kstr = sprintf(' sigma_{k=%d}', i);
    ktitles = cellfunc(@(s) [s, kstr], titles);
    for t = 1:T
        % TODO use http://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3-/content/munkres.m
        cequiv = oogreedymatching(gmms{t}.mu, isogmfit.mu); 
        assert(all(sort(cequiv) == 1:K));
        
        figure(1);
        subplot(T, K, sub2ind([K, T], i, t)); 
        imagesc(gmms{t}.sigma(:,:,cequiv(i))); title(ktitles{t});
        
        figure(2); hold on;
        subplot(1, K, sub2ind([K, 1], i, 1)); hold on;
        plot(gmms{t}.mu(cequiv(i), :)); 
    end
    subplot(1, K, sub2ind([K, 1], i, 1)); hold on;
    legend(titles);
    title(sprintf(' mu_{k=%d}', i)); 
end
figure(1); colormap gray;
figure(2); colormap gray;

%% resample from each gmm
T = numel(gmms);
for t = 1:T
    cequiv = oogreedymatching(gmms{t}.mu, isogmfit.mu);
    assert(all(sort(cequiv) == 1:K));
    samples{t} = gmms{t}.sample(10);
end

figure();
for t = 1:T
    subplot(1, T, t);
    imagesc(samples{t}); colormap gray;
end

    
%% Try some reconstruction
reconeigPatches = papago.recon(dswgmmm3, dspatches, W, 'eig', 5);
rX = {dspatches, reconeigPatches};
papago.visRecon(isopatches, W, rX, patchSize, 'titles', {'DS', 'EIG'});%, 'dorand', false);
err

%% Re-sample patches from distributions
% simulation using patches
% this way we'll have some patch-like-looking data
N = 10000;                   % number of patches to sample from GMM
D = size(isopatches, 2);    % dimensionality of the patch
K = 3;                      % number of clusters in GMM to learn from data

% learn and sample from a gmm from the nLearn random patches.
nlearn = min(size(isopatches, 1), max(10000, D)+1);
r = randperm(size(isopatches, 1));
Xlearn = isopatches(r(1:nlearn), :);
[X, mu, sigma, pi, cidx, sigmainv] = subspacetools.simgmm(N, D, K, 'frompatchesgmm', Xlearn);
W = max(maskpatches(1:N, :), 0.0001); % This is bad since some poitns will have all near-0 weights in 3x3
% W = unifrnd(0.0001, 1, [N, D]);
wgmmtrue = wgmm(X, W, K, mu, sigma, pi, sigmainv); 

% matlab's fitgmdist
gmmfit = fitgmdist(X, K, 'RegularizationValue', 0.00001, 'replicates', 5);
gmfit.sigma = gmmfit.Sigma;
gmfit.mu = gmmfit.mu;

% wgmm with W = 1
[fwgmmW1] = wgmmfit(X, W * 0 + 1, K, 'debug', false, 'sigmareg',  0.00001, 'replicates', 5);

% wgmm with W
[fwgmmW] = wgmmfit(X, W, K, 'debug', false, 'sigmareg',  0.00001, 'replicates', 5);

% wgmm with W with extra division
test1args = {'covarUpdateMethod', 'model3', 'muUpdateMethod', 'model3', 'logpUpdateMethod', 'model3'};
[fwgmmWm3] = wgmmfit(X, W, K, 'debug', false, test1args{:}, 'sigmareg',  0.0001, 'replicates', 20);
% Note: adding a large sigmareg screws up estimates.

% compare
gmms = {wgmmtrue, gmfit, fwgmmW1, fwgmmW, fwgmmWm3};
titles = {'true sigma', 'matlab gmm', 'W=1 wgmm m2', 'W wgmm m2', 'W wgmm m3'}; 
% gmms = {wgmmtrue, fwgmmWed};
% titles = {'true sigma', 'W wgmm modified'}; 

% show sigmas and mus
T = numel(gmms);
for i = 1:K; 
    kstr = sprintf(' sigma_{k=%d}', i);
    ktitles = cellfunc(@(s) [s, kstr], titles);
    for t = 1:T
        figure(1);
        subplot(T, K, sub2ind([K, T], i, t)); 
        
        [~, cequiv] = min(pdist2(gmms{t}.mu, wgmmtrue.mu));
        
        imagesc(gmms{t}.sigma(:,:,cequiv(i))); title(ktitles{t});
        
        figure(2); hold on;
        subplot(1, K, sub2ind([K, 1], i, 1)); hold on;
        plot(gmms{t}.mu(cequiv(i), :)); hold on;
    end
end
colormap gray;

