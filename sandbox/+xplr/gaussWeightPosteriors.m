%% xplr posterior assignments with weights
% extract iso patches from a location, and estimate a gaussian. Then assign posteriors.
% then, extract weight and ds patches, and assign posteriors via different methods
%   1. normal posterior of weights
%   2. using W.*Ds
%   3. normal posterior in rotated space.
% compare these results with iso posteriors

%% Load necessary data.

dsRate = 5;
usRate = 5;
K = 10;
location = LOC_VENTRICLE_EDGE;
locpad = ones(1,3)*1;
patchSize = ones(1,3)*9;
load([SYNTHESIS_DATA_PATH, 'adni_mfptrs_2015_11_08.mat']);

dsregmod = sprintf('brainDs%dUs%dReg', dsRate, usRate);
dsmaskregmod = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);
isoregmod = sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate);
[p, patches.layeridx, patches.volidx] = subspacetools.loadpatches('matfiles', mfptrs, ...
    {isoregmod, dsmaskregmod, dsregmod}, patchSize, location, locpad);

patches.iso = p{1};
patches.mask = p{2};
patches.ds = p{3};

%% estimate gaussians
wg = wgmmfit(patches.iso, patches.iso*0+1, K, 'replicates', 1,...
 'wgmmfields', struct('covarMergeMethod', 'none', 'sigmareg', 0.001));

%% reconstruct.
for k = 1:K
    % method 0
    isopost(:, k) = wg.logmvnpdf(patches.iso, wg.mu(k, :), wg.sigma(:,:,k));

    % method 1
    dspost(:, k) = wgmm.logmvnpdf(patches.ds, wg.mu(k, :), wg.sigma(:,:,k));

    % method 2
    wdspost(:, k) = wgmm.logmvnpdf(patches.ds .* patches.mask, bsxfun(@times, wg.mu(k, :), patches.mask), wg.sigma(:,:,k));

    % method 3
    warning('this needs a cleaner implementation of the conditional method, without for exampel the entire subject volume as input. and need posterior output.');
%     for reconSubj = 1:md.getNumSubjects()
%         subjdsnii = md.loadModality(sprintf('brainDs%dUs%d', dsRate, usRate), reconSubj);
%         subjdsvol = double(subjdsnii.img);
%         subjVolSize = size(subjdsvol);
%         subjDims = subjdsnii.hdr.dime.pixdim(2:4);
%         rSubjSpace = imref3d(subjVolSize, subjDims(2), subjDims(1), subjDims(3));
%         subjmask = nii2vol(md.loadModality(sprintf('brainDs%dUs%dMask', dsRate, usRate), reconSubj));
% 
%         % load subject data in atlas space (i.e. registered data)  
%         % necessary for gaussian conditional reconstruction
%         dsregnii = md.loadModality(sprintf('regBrainDs%dUs%d', dsRate, usRate), reconSubj);
%         atlVolSize = size(dsregnii.img); % size of volumes in atlas space
%         atlDims = dsregnii.hdr.dime.pixdim(2:4);
%         rAtlSpace = imref3d(atlVolSize, atlDims(2), atlDims(1), atlDims(3));
%         dsmaskregnii = md.loadModality(sprintf('regBrainDs%dUs%dMask', dsRate, usRate), reconSubj);
%         regtformmat = load(md.getModality(sprintf('brainDs%dUs%dRegMat', dsRate, dsRate), reconSubj));
% 
%         % get the locations in subject volume that correspond to voxel locations in the atlas volume 
%         locVolumeAtlas = getCorrespondingLoc(atlVolSize, regtformmat.tform, rAtlSpace, rSubjSpace, subjVolSize);
%         locVolumeSubject = getCorrespondingLoc(subjVolSize, regtformmat.tform.invert, rSubjSpace, rAtlSpace, atlVolSize);
% 
%         volRange = arrayfunc(@(x, p) (x: x + p - 1)', location, patchSize);
%         atlasPatch = dsregnii.img(volRange{:});
%         weightPatch = max(dsmaskregnii.img(volRange{:}), minW);
%         m = mean(atlasPatch(:));

%     end
    % rotpost
end

isopost = bsxfun(@minus, isopost, mean(isopost, 2));
dspost = bsxfun(@minus, dspost, mean(dspost, 2));
wdspost = bsxfun(@minus, wdspost, mean(wdspost, 2));

% compare methods
[~, m0] = max(isopost, [], 2);
[~, m1] = max(dspost, [], 2);
[~, m2] = max(wdspost, [], 2);

fprintf('simple ds posterior agreement with iso posterior: %1.3f\n', sum(m0 == m1) ./ numel(m0));
fprintf('w * ds posterior agreement with iso posterior: %1.3f\n', sum(m0 == m2) ./ numel(m0));
fprintf('rot(ds(w01)) posterior agreement with iso posterior: TODO\n');

%% Other visualizations
r32 = @(x) subspacetools.reshapeN3Dto2D(x, patchSize);
imagesc(r32(wg.mu));

figure();
for i = 1:5:size(isopost, 1); clf; plot(isopost(i, :)); hold on; plot(dspost(i, :)), hold on; plot(wdspost(i, :)); legend({'iso', 'ds', 'wds'}); pause; end


