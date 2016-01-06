%% setup
subspacesetup; 
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_VENTRICLE_EDGE; LOC_LEFT_CORTEX; 
locpad = o3 * 2;
maskthr = 0.75;

% error function
% TODO: compare different error functions
errfun = @(x, y) subspacetools.patcherror(x, y, @msd, 'subpatch', patchSize, o3*3, o3*7);
errfun = @(x, y) subspacetools.patcherror(x, y, @msd, 'gauss', patchSize, 1.5);

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);
nandspatches = dspatches;
nandspatches(maskpatches < maskthr) = nan;
nNans = sum(maskpatches(:) < maskthr);
fprintf('%d / %d are nan (%3.2f%%)\n', nNans, numel(maskpatches), nNans./numel(maskpatches));

%% analysis of linear interpolation patches 
% note these patches are extracted from whole volume interpolation, so they use more information
% than what's available in these patche columns.
errorLinInterp = errfun(dspatches, isopatches);

%% pca imputation
pcakrange = [2, 5, 10, 30]; [2:2:14, 20:5:30]; % TODO: just compute internally to catch 95% ?
nReps = 15;

st = statset('ppca'); st.Display = 'final'; st.MaxIter = 100;
opts = {'trueData', isopatches, 'realErrFun', errfun, 'maxIter', nReps};

pcaInitImputedPatches = cell(numel(pcakrange), 1);
pcaIsoImputedPatches = cell(numel(pcakrange), 1);
ppcaInitImputedPatches = cell(numel(pcakrange), 1);
pcaImputedPatches = cell(numel(pcakrange), 1);
ppcaImputedPatches = cell(numel(pcakrange), 1);

parfor ki = 1:numel(pcakrange)
    k = pcakrange(ki);
    pcaImputedPatches{ki} = pcaimpute(nandspatches, k, opts{:});
    pcaInitImputedPatches{ki} = pcaimpute(nandspatches, k, opts{:}, 'valinit', dspatches);
    
    % try with isopatches
    initpatches = isopatches; 
    nanisopatches = isopatches; nanisopatches(isnan(nandspatches)) = nan;
    pcaIsoImputedPatches{ki} = pcaimpute(nanisopatches, k, opts{:}, 'valinit', initpatches);
    
    % try ppca for missing values...
    [coeff,score,pcvar, mu, v, S] = ppca(nandspatches, k, 'Options', st);
    ppcaImputedPatches{ki} = bsxfun(@plus, score*coeff', mu); %S.Recon;
    
    % ppca with initialization. % TODO: try ppca_ implementation
    [~, ~, ~, ~, ~, S] = ppca(dspatches, k, 'Options', st); 
%     [~, ~, ~, ~, ~, S] = ppca(isopatches, k, 'Options', st); % try with iso patches to see how it works.
    [coeff, score, pcvar, mu, v, S] = ppca(nandspatches, k, 'Options', st, 'W0', S.W); 
    ppcaInitImputedPatches{ki} = bsxfun(@plus, score*coeff', mu); %S.Recon;
end
errorPCA = cellfun(@(x) errfun(x, isopatches), pcaImputedPatches);
errorIsoPCA = cellfun(@(x) errfun(x, isopatches), pcaIsoImputedPatches);
errorInitPCA = cellfun(@(x) errfun(x, isopatches), pcaInitImputedPatches);
errorPPCA = cellfun(@(x) errfun(x, isopatches), ppcaImputedPatches);
errorInitPPCA = cellfun(@(x) errfun(x, isopatches), ppcaInitImputedPatches);



%% gmm imputation
% TODO: idea: since gmm is slow for real data with needed large C, maybe do pca first, then gmm on
% that, then go back to full space ?
gmmpcakrange = [2, 5, 10, 30]; % TODO: just compute internally to catch 95% ?
crange = [3, 20, 50, 100]; %[2, 3, 5, 10, 20, 100];
nReps = 15;
opts = {'trueData', isopatches, 'realErrFun', errfun, 'maxIter', nReps};
niifiles = arrayfunc(@(x) sprintf('subspace/data/bucknerNiis_5_%d.mat', x), 2:5);

% HGMM options
dsRate = 5;
usRates = 2:dsRate;
patchFactors = 2*round(linspace(3, 4, numel(usRates))) + 1;

% gmm methods (k=0 means no PCA)
gmmImputedPatches = cell(numel(crange), numel(gmmpcakrange));
gmd = cell(numel(crange), numel(gmmpcakrange));
hgmmImputedPatches = cell(numel(crange), numel(gmmpcakrange));
hgmd = cell(numel(crange), numel(gmmpcakrange));
wgmmImputedPatches = cell(numel(crange), numel(gmmpcakrange));
wgmd = cell(numel(crange), numel(gmmpcakrange));
for ci = 1:numel(crange)
    c = crange(ci);
    parfor ki = 1:numel(gmmpcakrange)
        k = gmmpcakrange(ki);
        warning off backtrace; 
        
        % hierarchical GMM
        [vols, hgmmImputedPatches{ci, ki}, hgmd{ci, ki}] = ...
            hgmmimpute(niifiles, dsRate, usRates, patchFactors, location, locpad, c, niis.iso, k, nReps);
        
        % TODO: would be great to avoid volinit! volinit is misleading :( 
        % maybe init with ppca results?
        [gmmImputedPatches{ci, ki}, gmd{ci, ki}] = ...
            gmmimpute(nandspatches, c, k, 'valinit', dspatches, opts{:});
        
        % try the weight.
        [wgmmImputedPatches{ci, ki}, wgmd{ci, ki}] = ...
            gmmimpute(dspatches, c, k, 'valinit', dspatches, 'valwts', maskpatches, opts{:});
        
        % ToTry: initiate with S.Recon (do ppca --> S.recon outside of clusters?)
    end
end
errorGMMPCA = cellfun(@(x) errfun(x, isopatches), gmmImputedPatches);
errorWGMMPCA = cellfun(@(x) errfun(x, isopatches), wgmmImputedPatches);
errorHGMMPCA = cellfun(@(x) errfun(x, isopatches), hgmmImputedPatches);
        

% visualize result.
i = 500;
t = @(x) subspacetools.reshape3Dto2D(x(i, :), patchSize);
tpatches = cellfunc(t, {isopatches, maskpatches, dspatches, abs(dspatches-isopatches), gmmImputedPatches{ci, ki}, abs(gmmImputedPatches{ci, ki}-isopatches), wgmmImputedPatches{ci, ki}, hgmmImputedPatches{ci, ki}});
%view2D(tpatches); colormap gray;
figure(); imagesc(cat(1, tpatches{:})); colormap gray;

%% visualization of results
figuresc(); hold on;

% display linear interpolation4
plot([min([pcakrange, gmmpcakrange]), max([pcakrange, gmmpcakrange])], errorLinInterp * [1, 1]); 

% display PCA
% plot(pcakrange, errorPCA, '*-'); hold on; 
% plot(pcakrange, errorInitPCA, 'x-'); hold on; 
% plot(pcakrange, errorPPCA, '+-'); hold on; 
% plot(pcakrange, errorInitPPCA, 'd-'); hold on; 
plot(gmmpcakrange, errorGMMPCA, 'o-'); hold on; 
plot(gmmpcakrange, errorWGMMPCA, 'd-'); hold on; 
plot(gmmpcakrange, errorHGMMPCA, 's-'); hold on; 
xlabel('K'); ylabel('mean msd error');
gmmpcaTitles = arrayfunc(@(x) sprintf('GMM-PCA imput %d clust', x), crange);
wgmmpcaTitles = arrayfunc(@(x) sprintf('WGMM-PCA imput %d clust', x), crange);
hgmmpcaTitles = arrayfunc(@(x) sprintf('HWGMM-PCA imput %d clust', x), crange);
pcaTitles = {'PCA', 'PCA w/ lin init', 'PPCA', 'PPCA w/ lin init'};
legend([{'(vol) linear interp'} ...
    gmmpcaTitles, wgmmpcaTitles, hgmmpcaTitles]);
title('inputation results');
