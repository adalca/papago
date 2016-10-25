% Simulate data and test aspects of model4-wgmm

% initialize
setup

%% parameters
nSimSamples = 4750;

atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
atlPatchSize = ones(1, 3) * 6; 
atlLoc = [69, 51, 43]; % scale 3

atlPatchSize = ones(1, 3) * 5;
atlLoc = [20, 25, 35];
atlLoc = [20, 25, 35];

gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 7; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
traindataset = 'adni';
testdataset = 'adni';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad);

% modalities
ds = 5;
us = 2;
isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
dsSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Reg', ds, us);
dsSubjInAtlMaskMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2MaskReg', ds, us);
dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);

dsSubjMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2', ds, us);
dsSubjMaskMod = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Mask', ds, us);
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);


% isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
% dsSubjInAtlMod = sprintf('brainDs%dUs%dReg', ds, us);
% dsSubjInAtlMaskMod = sprintf('brainDs%dUs%dRegMask', ds, us);
% dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
% dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);
% 
% dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
% dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);
% isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load various data columns from training set
fnames = fullfile(SYNTHESIS_DATA_PATH, traindataset, 'md', [sys.usrname, '_restor_md_*']);
trainmd = loadmd(fnames);
[bucknerIsoPatchCol, layeridx, volidx] = ...
    subspacetools.md2patchcol(trainmd, isoSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsMaskPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjInAtlMaskMod, atlPatchSize, atlLoc, patchColPad);


% load selected ADNI subject volumes
fnames = fullfile(SYNTHESIS_DATA_PATH, testdataset, 'md', [sys.usrname, '_restor_md_*']);
testmd = loadmd(fnames);
dsSubjNii = testmd.loadModality(dsSubjMod, reconSubj);
dsSubjVol = double(dsSubjNii.img);
dsSubjWeightVol = logical(testmd.loadVolume(dsSubjMaskMod, reconSubj));
dsSubjInAtlNii = testmd.loadModality(dsSubjInAtlMod, reconSubj);
dsSubjInAtlMaskVol = testmd.loadVolume(dsSubjInAtlMaskMod, reconSubj);
subjInAtlTform = load(testmd.getModality(dsSubjInAtlMatMod, reconSubj));

% prepare necessary inputs for conditional-based reconstruction
subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
tform = subjInAtlTform.tform;
atlVolSize = size(dsSubjInAtlNii.img);
subjLoc2AtlSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims);
atlLoc2SubjSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims, 'backward');
extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);

%% verify data
err = abs(bucknerIsoPatchCol - bucknerDsPatchCol);
plot(bucknerDsMaskPatchCol(:), err(:), '.');

l = linspace(0, 1, 100);
z = kernelRegress(bucknerDsMaskPatchCol(:), err(:), l, 0.05);
hold on; plot(l, z, '.-');

axislabels('patch weights', 'abs error', '');
legend({'data', 'kernel regression'});

%% new D function using adjacency
q = size2ndgridvec(atlPatchSize);
Adj = exp(-pdist2(q, q));
fn1 = @(w) diag((-log(max(w, eps))).^2);
fn2 = @(w) Adj * diag((-log(w)).^2);
fn3 = @(w) diag((-log(w))) * Adj * diag((-log(w))) * 0.06;% original factor was 0.0006
fn4 = @(w) (-log(w))' * (-log(w)) .* Adj * 0.06;
% view2D({D1, D2; D3 D4}, 'titles', {'D=diag(-log(W)).^2', 'Adj*D', 'sqrt(D)*Adj*sqrt(D)', '(W''W)*Adj'}); colormap gray;
fnchoice = fn1;

%% learn a gmm
gmmopt = statset('Display', 'iter', 'MaxIter', 20, 'TolFun', 0.001);

% compute the gaussian mixture model. TODO: should use wgmmfit with model0
tic
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
gmdistIso = fitgmdist(X, gmmK, 'regularizationValue', regVal, 'replicates', 3, 'Options', gmmopt);
fprintf('gmm took %3.3f sec\n', toc);
gmmIso = wgmm.gmdist2wgmm(gmdistIso);

isopost = gmdistIso.posterior(X);
[~, mi] = max(isopost, [], 2);

%% simulate destorying some of the data according to the model
gmmIsoK1 = wgmm(gmmIso.mu(1, :) * 0, gmmIso.sigma(:, :, 1), 1);
% tmp = rand(size(Xt, 2),size(Xt, 2)); 
% sigma = tmp'*tmp; 
% gmmIsoK1 = wgmm(gmmIso.mu(1, :) * 0, sigma, 1);
[Xt, cl] = gmmIsoK1.sample(nSimSamples);
W = rand(size(Xt));
W(W < 0.1) = 0.1;
X0 = zeros(size(Xt));
for i = 1:size(Xt, 1)
    w = W(i, :);
    Di = fnchoice(w);
    X0(i, :) = mvnrnd(Xt(i, :), Di);
end


%% "true" data

k = 1;

Xt = bucknerIsoPatchCol(mi == k, :);
Xtmean = mean(Xt, 2);
Xt = bsxfun(@minus, Xt, Xtmean);
W = bucknerDsMaskPatchCol(mi == k, :);
W(W < 0.00001) = 0.00001;

X0 = bucknerDsPatchCol(mi == k, :);
X0 = bsxfun(@minus, X0, Xtmean);
x0err = msd(Xt, X0, 2);

nSimSamples = size(Xt, 1);
gmmIsoK1 = wgmm(gmmIso.mu(k, :), gmmIso.sigma(:, :, k), 1);

%% try several sigmas

% get all poitns which belong to cluster 1
truesigma = gmmIsoK1.sigma;
truemu = gmmIsoK1.mu;

% recover from true samples Xt, using truemu
numer = 0; denom = 0;
for i = 1:nSimSamples
    x = Xt(i, :);
    df = (x - truemu);
    dfd = df(:) * (df(:))';
    
    numer = numer + dfd;
    denom = denom + 1;
end
sigmarecon{1} = denom \ numer; 
sigmatitles{1} = 'Xt';

% recover from x0 samples, true mean
numer = 0; denom = 0;
for i = 1:nSimSamples
    x = X0(i, :);
    df = (x - truemu);
    dfd = df(:) * (df(:))';
    
    numer = numer + dfd;
    denom = denom + 1;
end
sigmarecon{2} = denom \ numer; 
sigmatitles{2} = 'x0';

% "recover" using hack 1
itersigma = sigmarecon{2};
for j = 1:1
    numer = 0; denom = 0;
    for i = 1:nSimSamples
        w = W(i, :); 
        x = X0(i, :);
        df = (x - truemu);
        dfd = df(:) * (df(:))';

        Di = fnchoice(w);
        sg = itersigma + Di;

        numer = numer + (sg \ dfd);
        denom = denom + inv(sg);
    end
    itersigma = denom \ numer;
end
sigmarecon{3} = itersigma; 
sigmatitles{3} = 'hack1';

% "recover" using hack 2
numer = 0; denom = 0;
for i = 1:nSimSamples
    w = W(i, :); x = X0(i, :);
    w(w == 1 ) = 0.99999;
    df = (x - truemu);
    dfd = df(:) * (df(:))';

    Di = fnchoice(w);

    numer = numer + (dfd / sg);
    denom = denom + inv(sg);
end
sigmarecon{4} = numer / denom; 
sigmatitles{4} = 'hack2';

% recover via hack 3 (the cornell hack)
numer = 0; denom = 0;
for i = 1:nSimSamples
    w = W(i, :); 
    x = X0(i, :);
    w(w == 1 ) = 0.99999;
    df = (x - truemu);
    dfd = df(:) * (df(:))';

    Di = fnchoice(w);

    numer = numer + dfd  - Di;
    denom = denom + 1;
end
sigmarecon{5} = denom \ numer; 
sigmatitles{5} = 'hack3 (wornell)';

% zfn = @(s) model5der(s, X(cl1idx, :), W(cl1idx, :), truemu, ones(nSimSamples, 1), fnchoice);
% sigmav3 = fzero(zfn, truesigma);


% "recover" using hack 1 with subtracting meanD
itersigma = sigmarecon{2};
for j = 1:1
    numer = 0; denom = 0; meanD = 0;
    for i = 1:nSimSamples
        w = W(i, :); x = X0(i, :);
        df = (x - truemu);
        dfd = df(:) * (df(:))';

        Di = fnchoice(w);
        sg = itersigma + Di;

        numer = numer + (sg \ dfd);
        denom = denom + inv(sg);
        
        meanD = meanD + Di;
    end
    meanD = meanD ./ nSimSamples;
    itersigma = (denom \ numer);
end 
sigmarecon{6} = itersigma - meanD; 
sigmatitles{6} = 'hack1v2';

sigmarecon{7} = sigmarecon{3} - 0.25 * meanD; 
sigmatitles{7} = 'hack1v3';

% "recover" using hack 8
itersigma = sigmarecon{2};
for j = 1:1
    numer = 0; denom = 0; meanD = 0;
    for i = 1:nSimSamples
        w = W(i, :); x = X0(i, :);
        df = (x - truemu);
        dfd = df(:) * (df(:))';

        Di = fnchoice(w);
        sg = itersigma + Di;

        numer = numer + (sg \ dfd);
        denom = denom + inv(sg) + inv(sg) * dfd * inv(itersigma) * Di * inv(itersigma);
        
        meanD = meanD + Di;
    end
    meanD = meanD ./ nSimSamples;
    itersigma = (denom \ numer);
end 
sigmarecon{8} = itersigma - meanD/8; 
sigmatitles{8} = 'hack8';

% "recover" using hack 9
itersigma = sigmarecon{2};
for j = 1:1
    numer = 0; denom = 0; meanD = 0;
    for i = 1:nSimSamples
        w = W(i, :); x = X0(i, :);
        df = (x - truemu);
        dfd = df(:) * (df(:))';

        Di = fnchoice(w);
        sg = itersigma + Di;
        

        numer = numer + (sg \ dfd);
        t0 = Di / itersigma;
        t1 = inv(itersigma) * t0;
        t2 = t1 * t0;
        t3 = t2 * t0;
        denom = denom + inv(sg) + inv(sg) * dfd * t1 - inv(sg) * dfd * t2; % + inv(sg) * dfd * t3;
        
        meanD = meanD + Di;
    end
    meanD = meanD ./ nSimSamples;
    itersigma = (denom \ numer);
end 
sigmarecon{9} = itersigma; 
sigmatitles{9} = 'hack9';

%% visualize
vs = volstats(truesigma);
view2D([truesigma, sigmarecon], 'titles', ['orig', sigmatitles], 'caxis', [vs.min, vs.max]); colormap gray;
% 
% figure(); imagesc([truesigma, sigmarecon{[3, 5]}]);

%% reconstruction of X
wienerfilt = @(s,Di,y,mu) (s / (s + Di)) * (y(:) - mu(:)) + mu(:);

xrecon = {};
for i = 1:size(X0, 1)
    y = X0(i, :);
    w = W(i, :);
    Di = fnchoice(w);
    
    for si = 1:numel(sigmarecon)
        s = sigmarecon{si};
        xrecon{si}(i, :) = wienerfilt(s, Di, y, truemu);
    end
    
    truexrecon(i, :) = wienerfilt(truesigma, Di, y, truemu);
end

%%
msdq = cellfunc(@(x) msd(x, Xt, 2), xrecon);
madq = cellfunc(@(x) mean(abs(x - Xt), 2), xrecon);
fprintf('%20s %20s %20s\n', 'sigma-type','mean-msd', 'median-msd');
for i = 1:numel(sigmarecon)
    fprintf('%20s %20.7f %20.7f %20.7f\n', sigmatitles{i}, mean(msdq{i}), median(msdq{i}), median(madq{i}));
end
msdq = msd(truexrecon, Xt, 2);
madq = mean(abs(truexrecon - Xt), 2);
fprintf('%20s %20.7f %20.7f %20.7f\n', '"truesigma"', mean(msdq), median(msdq), median(madq));
msdq = msd(X0, Xt, 2);
madq = mean(abs(X0 - Xt), 2);
fprintf('%20s %20.7f %20.7f %20.7f\n', 'actual-x0', mean(msdq), median(msdq), median(madq));


%% reconstruct subvolume
[quiltedSubvolSio, reconLoc, cntvol] = papago.subvolRecon(gmmIso, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


[quiltedSubvolM5, reconLoc, cntvol] = papago.subvolRecon(wg, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);









%%
nsSamples = 10;
r = randsample(size(X0, 1), nsSamples);

z = {};
for i = r(:)'
    
    x0_ = X0(i, :);
    xt_ = Xt(i, :);
    stt_ = xt_(:) * (xt_(:))';
    s00_ = x0_(:) * (x0_(:))';
    w_ = W(i, :);
    z = [z, abs(stt_ - s00_), w_(:) * w_(:)'];
end

view2D(z)

% 
%% E-step (mixture of wiener filters) check

% try just making a new wgmm, and do estep with the correct sigma, mu
wg = wgmm(gmmIso.mu, gmmIso.sigma, gmmIso.pi); % ! with defaults set to model 4.
[~, ~, Xhat] = wg.estep(X0, W);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));

%% run full model4-wgmm on this simulated data.
% estimate gmm.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3); % ! with defaults set to model 4.

% estimate Xhat and cluster assignments
[~, gammank, Xhat] = wg.estep(X0, W);
[~, cidxhat] = max(gammank, [], 2);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));

figure();
hist(cidxhat, 1:gmmK);

% TODO: use wgmm model 3 to get sigma, mu, and weiner filter to get xhat, and compare

%% simulate data from cluster, but with real weights
[Xt, cl] = gmmIso.sample(nSimSamples);

W = bucknerDsMaskPatchCol(randsample(size(bucknerDsMaskPatchCol, 1), nSimSamples), :);
W = max(W, 0.000001);

X0 = zeros(size(Xt));
for i = 1:size(Xt, 1)
    w = W(i, :);
    Di = diag((-log(w)).^2);
    X0(i, :) = mvnrnd(Xt(i, :), Di);
end

% estimate gmm.
wg = wgmm.fit(X0, W, gmmK, 'replicates', 3); % ! with defaults set to model 4.

% estimate Xhat and cluster assignments
[~, gammank, Xhat] = wg.estep(X0, W);
[~, cidxhat] = max(gammank, [], 2);
xhaterr = msd(Xt(:), Xhat(:));

figure(); % xhat vs x0 error
subplot(221); plot(Xt, X0, 'b.'); title(sprintf('original err, msd: %3.2f', x0err));
subplot(222); plot(Xt, Xhat, 'b.'); title(sprintf('xhat err, msd: %3.2f', xhaterr));
subplot(223); histogram(msd(Xt, X0, 2)); title(sprintf('original err hist'));
subplot(224); histogram(msd(Xt, Xhat, 2)); title(sprintf('xhat err hist'));


