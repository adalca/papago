% explore the weight to noise relationship resulting from downsampled data.

% initialize
setup

%% parameters
nSimSamples = 1000;

atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
% atlPatchSize = ones(1, 3) * 5; 
% atlLoc = [20, 25, 35];
gmmK = 15; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 3; %1, 3
patchColPad = ones(1, 3) * 1;

% train and test datasets
traindataset = 'adni';
testdataset = 'buckner';

% recon params
keepr = 1;
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

subvolSize = atlPatchSize + 5;

% modalities
ds = 5;
us = 5;
isoSubjInAtlMod = sprintf('brainIso2Ds%dUs%dsizeReg', ds, us);
dsSubjInAtlMod = sprintf('brainDs%dUs%dReg', ds, us);
dsSubjInAtlMaskMod = sprintf('brainDs%dUs%dRegMask', ds, us);
dsSubjInAtlMatMod = sprintf('brainDs%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
dsInterpSubjInAtlMod = sprintf('brainDs%dUs%dInterpReg', ds, us);

dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);


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


%% analyze
err = abs(bucknerIsoPatchCol - bucknerDsPatchCol);
plot(bucknerDsMaskPatchCol(:), err(:), '.');

l = linspace(0, 1, 100);
z = kernelRegress(bucknerDsMaskPatchCol(:), err(:), l, 0.05);
hold on; plot(l, z, '.-');

%% sample from the data, and compute weight covariance, error covariance.
nSamples = 1000;

r = randsample(size(bucknerIsoPatchCol, 1), nSamples); 

dall = [];
scall = [];
for i = 1:nSamples
    % extract sampled patches
    w = bucknerDsMaskPatchCol(r(i), :);
    Xi = bucknerIsoPatchCol(r(i), :);
    X0 = bucknerDsPatchCol(r(i), :);
    x = (X0 - Xi);

    % weight covariance. square matrix.
    D = w'*w; 
    
    % error covariance. square matrix.
    sc = x' * x;
    
    % add to long vector. this will be slow since I'm not pre-allocating.
    dall = [dall; D(:)];
    scall = [scall; sc(:)];
end

%% compute means and standard deviations per weight-bin 
nBins = 100;
l = linspace(eps, 1-eps, nBins);
[~, ~, assn] = histcounts(dall, l); 

for i = 1:nBins
    data = scall(assn == i);
    means(i) = nanmean(data);
    stds(i) = nanstd(data);
    i
end    
    
%% visualize results.
rShow = randsample(size(dall, 1), 1000000);
figure(); 
subplot(1,2,1);
plot(dall(rShow), scall(rShow), '.'); hold on;
axislabels('ww''(:)', '(x0-xt)(x0-xt)''(:)', 'ds5us5 random sample of patches analysis');
plot(linspace(0, 1, nBins), stds, '.-k');
plot(linspace(0, 1, nBins), 2.5 .* stds, '.-r');
plot(l, -log(l)*0.0006, '-y');
xlim([0,1])
ylim(0.02*[-1,1])
legend({'data', 'std', '2.5std', '-log(w^2)*0.0006'})

subplot(1,2,2);
plot(linspace(0, 1, nBins), stds, '.-k');
data = scall(assn == 10);
hist(data, 1000); title(sprintf('w^2 = %3.2f', l(10))); xlim([-1,1]*0.005);
