% add path
subspacesetup;

% load data
load('data/bucknerNiis_5_5.mat'); 

% parameters
patchSize = ones(1, 3) * 5;
patchSizeIn2D = [patchSize(1), patchSize(2)*patchSize(3)];
location = [56, 72, 92];
locpad = ones(1, 3) * 2;
nClust = 5; 


% load in all the volumes into a cell array
nVols = length(niis.iso);
for i=1:nVols; 
    vols{i} = niis.iso{i}.img;  
end

% extract subvolumes
volRange = arrayfunc(@(x, d, p) (x - d: x + d + p - 1)', location, locpad, patchSize);
vols = cellfunc(@(ns) ns(volRange{:}), vols);
    
%% ========== REGULAR GMM ============= %%

% prepare data (patches)
libs = patchlib.vol2lib(vols, patchSize);
patches = cat(1, libs{:});

gmd = gmdistribution.fit(patches, nClust, 'Replicates', 5, 'Regularize', 0.00001); 

% sample from regular gmm
samplesRegular = gmmSampleGrid(gmd, patchSizeIn2D, [10, 2]);

%% ========== ONLINE GMM ============= %%

filename = 'onlineGMMPatchPrior2'; 
MiniBatchSize = round(0.1*size(patches,1));
nIter = 1000; 

%% learn model from training data
NewGMM = OnlineGMMEM(nClust,@(N) removeDC(Rand3DPatchesFromImagesCell(N,patchSize,vols)),nIter,MiniBatchSize,filename,2,0.9,MiniBatchSize,true);

% sort output
[NewGMM.mixweights,inds] = sort(NewGMM.mixweights,'descend');
NewGMM.covs = NewGMM.covs(:,:,inds);

gmdOnline.Sigma = NewGMM.covs; 
gmdOnline.mu = zeros(nClust, prod(patchSize)); 

% sample from online gmm
samplesOnline = gmmSampleGrid(gmdOnline, patchSizeIn2D, [10, 2]);

