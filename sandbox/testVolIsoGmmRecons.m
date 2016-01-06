%% Description 
% test the volume-wise (i.e. joinging patch-wise) reconstruction methods using the gmm parameters 
% from running matlab's standard GMM method on the iso data. 
import subspacetools.loadpatches;

%% General setup
o3 = ones(1, 3);
locpad = o3 * 1;

% data parameters
dsRate = 5; % loading
opts.load.atldstype = 'ds'; % low resolution data to run on: 'bisow' or 'ds'
opts.load.trainmfptrsfile = [SYNTHESIS_DATA_PATH, 'adni_mfptrs_2015_11_08.mat'];
opts.load.testmfptrsfile = [SYNTHESIS_DATA_PATH, 'buckner_mfptrs_2015_11_08.mat'];
opts.load.trainmdfile = [SYNTHESIS_DATA_PATH, 'adni_restor_md_2015_11_08.mat'];
opts.load.testmdfile = [SYNTHESIS_DATA_PATH, 'buckner_restor_md_2015_11_08.mat'];
opts.load.method = 'matfiles'; % 'volumes' or 'matfiles' is most likely
% opts.load.method = 'volumes';

% gmm and reconstruction parameters
gmmMethods = {'iso-gmm', 'iso-gmm_exclsubj', 'iso-wgmm_exclsubj'}; %, 'ds-gmm', 'ds-wgmm'};
gmmMethods = {'iso-gmm', 'ds-wgmm'};
K = 10; % used to vary this by location, but it's a bit difficult now with the slower volume loading
weightfact = prod(locpad*2+1) * 15; % sigma-recon weight threshold
regval = 1e-5;
minW = 1e-9;
gmmopt = statset('Display', 'final', 'MaxIter', 5, 'TolFun', 0.01);
gmmargs = {'RegularizationValue', regval, 'replicates', 3, 'Options', gmmopt};
wgmmargs = {'sigmareg', regval, 'replicates', 3, 'wgmmfields', struct('sigmaopt', {{weightfact}}, 'sigmareg', regval), 'TolFun', 0.01};

% reconstruction parameters
reconSubj = 4;
reconSigmaMethod = 'inverse'; %'forward', 'inverse'
reconRegVal = 1e-4; % regulaization to the diagonal of newSigma

pgmm = {};

%% scale-specific parameters
doRates = 2:dsRate;
patchSizes = arrayfunc(@(x) o3 * x, [5, 7, 9, 9]);
patchOverlaps = arrayfunc(@(x) o3 * x, [3, 3, 5, 5]);
cropranges = {[0.25, 0.75], [0.3, 0.7], [0.3, 0.5], [0.25, 0.55]};

%% do the core work
usirange = 1:numel(doRates);
usirange = numel(doRates);
for usi = usirange
    usRate = doRates(usi);
    patchSize = patchSizes{usi};
    patchOverlap = patchOverlaps{usi};
    
    % TODO: need to fix this.
    atlVolSize = size(nii2vol(loadNii(eval(sprintf('BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate)))));
    reconVolRange{1} = cropranges{usi}(1) * atlVolSize;
    reconVolRange{2} = cropranges{usi}(2) * atlVolSize; % subvolume to reconstruct. 
    croprange = arrayfunc(@(x, y) x:y, o3, reconVolRange{2} + patchSize + locpad);
    
    if usi == usirange(1)
        mods.dsreg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
    else
        mods.dsreg = sprintf('brainDs%dUs%dReg_rebuilt', dsRate, usRate);
    end
        
    mods.dsmaskreg = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);
    mods.isoreg = sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate);
    % bisoregmod = 'brigidRegBrain5_sigma2';
    
    % TODOs: 
    % 1. fix sigmamerge fact passing. Done (?)
    % 2. fix md to have _rebuilt modalities for ds_us
    % 3. save gmms, and load them as needed, building from previous gmm as init.
    %   have option for this. 'learn with init' or something like that.
    %
    % N. try global gmm without mean
    % N. try iso-gmm on ADNI and do on buckner
    % N. a. re-process buckner.
    %
    % do all recon methods, and always test them patch-by-patch and for whole volume!!!
    
    testVolIsoGmmRecons_Core;
end

%% combine patches into volume
% the problem is we can't use patchlib.stackPatches since we don't have the patches necessarily
% landing on a grid, due to rounding errors (I haven't checked this, just assuming...)
clear wmeanVol meanVol wVol;
ip = @subspacetools.irregularPatches2vol;
for t = 1:numel(gmmMethods)
    % compute weights based on how far they are from the nan margins.
    patchweights = cellfunc(@(x) double(bwdist(isnan(x))), reconPatches1{t});

    wmeanVol{t} = ip(reconSubs1{t}, reconPatches1{t}, 'volSize', size(subjdsvol));
    meanVol{t} = ip(reconSubs1{t}, reconPatches1{t}, 'volSize', size(subjdsvol), 'weightPatches', patchweights);
    wVol{t} = ip(reconSubs1{t}, patchweights, 'volSize', size(subjdsvol));
    
%     q = dimsplit(1, subgridlocs); q = q(1:numel(weightrecon1{t}));
%     q(end) = []; wr = weightrecon1{t}; wr(end) = [];
%     wr = cellfunc(@(x) reshape(x, patchSize), wr);
%     weightreconVol{t} = ip(q', wr, 'volSize', atlVolSize);
end

% check that in non-nan areas, the regenerated volumes match the original volume where mask is 1.
for t = 1:numel(gmmMethods)
    v = wmeanVol{t};
    maps = ~isnan(v) & subjmask;
    
    assert(all(isclose(wmeanVol{t}(maps), subjdsvol(maps))));
    assert(all(isclose(meanVol{t}(maps), subjdsvol(maps))));
end

%% visualization 3D
t = 1;
% another way to buld the volume, just average to compare//check
filledvol = nan(size(subjdsnii.img));
for i = find(~cellfun(@isempty, reconPatches1{t}))
    thisvol = setSubvolume(filledvol*nan, reconPatches1{t}{i}, reconSubs1{t}{i}, reconSubs1{t}{i} + size(reconPatches1{t}{i}) - 1);
    filledvol = nanmean(cat(4, thisvol, filledvol), 4);
end

filledmask = (~isnan(filledvol))*1;
% view3Dopt(meanVol, wmeanVol, wVol, filledvol, double(dsnii.img .* filledmask), double(isonii.img .* filledmask))

%% visualize in 2D (middle of compputed volume)
[croppedVol, cropMask, cropArray, bBox] = boundingBox(filledmask);

subjisonii = testmd.md.loadModality(sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate), reconSubj);
assert(all(subjVolSize == size(subjisonii.img)));

% isoregnii = loadNii(sprintf('data/%s/%s_brain_reg_adni61_brain.nii.gz', bpref, bpref)); % ds
% assert(all(size(isoregnii.img) == atlVolSize));

figure();
% get centroid, range
centroid = round(cellfun(@median, cropArray));
sliceRange = cropArray;
sliceRange{1} = centroid(1);
% sliceRange{2} = cropArray{2}(15)
% want to show: isonii, dsmasknii, dsnii \\ recons 1..t
T = numel(gmmMethods);
vf = @(x) squeeze(x(sliceRange{:}) .* filledmask(sliceRange{:}));
m = max(T, 3);
subplot(2, m, 1); imagesc(vf(subjisonii.img), [0, 0.75]); axis off; colormap gray; title('iso');
subplot(2, m, 2); imagesc(vf(subjmask), [0, 0.75]); axis off; colormap gray; title('dsmask');
subplot(2, m, 3); imagesc(vf(subjdsnii.img), [0, 0.75]); axis off; colormap gray; title('ds');
%subplot(2, T, 4); imagesc(vf(blurisowvols{subject}), [0, 0.75]); axis off; colormap gray; title('bisow');

for t = 1:T
    subplot(2, m, m + t); imagesc(vf(meanVol{t}), [0, 0.75]); axis off; colormap gray; title(gmmMethods{t});
end

%% visualization in 3d
subjdsregion = double(subjdsnii.img .* filledmask);
subjmaskregion = double(subjmask .* filledmask);
subjisoregion = double(subjisonii.img .* filledmask);

vols = {subjisoregion, subjmaskregion, subjdsregion, wmeanVol{:}};

crange = arrayfunc(@(x, y) x:y, reconSubs1{1}{1}, reconSubs1{1}{end-1} + patchSize + locpad);
cvols = cellfunc(@(x) x(crange{:}), vols);
cvols = cellfunc(@(x) within([0, 0.4], x), cvols);
view3Dopt(cvols{:}, cvols{end}*0);


