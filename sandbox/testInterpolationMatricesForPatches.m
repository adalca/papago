% params
atlLoc = [100, 100, 100];
atlPatchSize = [15, 15, 15];

% path setup
subjpath = 'D:/Dropbox (MIT)/Research/patchSynthesis/data/adni/tmp/101411/';

% load in subject_space ds5us5, atlas_space ds5us5, R mat file.
% choose random patch location and size
subj = nii2vol([subjpath, '101411_ds5_us5.nii.gz']);
atl = nii2vol([subjpath, '101411_ds5_us5_reg.nii.gz']);
% loads 'atl2subjR' 'subj2atlR' 'subjLoc2AtlSpace' 'atlLoc2SubjSpace'
load([subjpath, '101411_cor_2_ds5_us5_size.mat']);

% extract patch from atlas space --> pA
pA = cropVolume(atl, atlLoc, atlLoc + atlPatchSize - 1);

% get location of same patch in subject space. 
% extract patch from subject space --> pS
[~, subjPatchMins, subjPatchSize] = paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);
pS = cropVolume(subj, subjPatchMins, subjPatchMins + subjPatchSize - 1);

% get R (interpolation A->S) matrix for this patch
% rotate atlas space patch to subject space -->pAS
[RAS, atlMaskAS, subjMaskAS, atlRCropInd, subjRCropInd, subjWeight] = vol2subvolInterpmat(atl2subjR, atlLoc2SubjSpace, size(subj), atlLoc, atlPatchSize);
pAS = reshape(RAS * pA(:), size(pS));
pAS(~subjMaskAS) = nan; % note: y_i^S should be just the entries inside the mask.

% get /Gamma (interpolation S->A) matrix for this patch
% rotate subj patch to atl patch --> PSA
RSA = subj2atlR(atlRCropInd, subjRCropInd);
[~, ~, ~, ~, ind, atlWeight] = vol2subvolInterpmat(subj2atlR, subjLoc2AtlSpace, size(atl), subjPatchMins, subjPatchSize);
pSAtmp = reshape(RSA * pAS(:), size(pA));
pWeight1 = reshape(RSA * subjWeight(:), size(pA));

[~, ia, ~] = intersect(ind, atlRCropInd);
atlWeightCropped = reshape(atlWeight(ia), size(pA));

pWeight2 = reshape(RSA * subjWeight(:), size(pA)) .* atlWeightCropped;




% visually compare PA with PSA, and PS with PAS.






