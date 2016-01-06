%% parameters

subspacesetup;

patchSize = [9 9 9];
location = LOC_VENTRICLE_EDGE-5; LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;

regVal = 1e-4; %regulaization to the diagonal of newSigma

gmmmethod = 'groundtruth'; %'compute', 'load', 'groundtruth'
sigmaMethod = 'inverse'; %'forward', 'inverse'
subject = 4; %1, 3

%% load volumes
% we'll use some local, small subset, of the buckner data for now.

% TODO: should move this to using the existing md structure, 
% once both katie and adrian have access to the same large data.
bpref = sprintf('buckner%02d', subject);

% subject space
isonii = loadNii(sprintf('data/%s/%s_brain_cropped5.nii.gz', bpref, bpref)); % iso
dsnii = loadNii(sprintf('data/%s/%s_brain_downsampled5_reinterpolated5.nii.gz', bpref, bpref)); % ds 
dsmasknii = loadNii(sprintf('data/%s/%s_brain_downsampled5_reinterpolated5_dsmask.nii.gz', bpref, bpref));

% atlas space (i.e. registered data)
dsregnii = loadNii(sprintf('data/%s/%s_brain_downsampled5_reinterpolated5_reg.nii.gz', bpref, bpref)); % ds
load(sprintf('data/%s/%s_brain_downsampled5_reinterpolated5_reg.mat', bpref, bpref));
dsmaskregnii = loadNii(sprintf('data/%s/%s_brain_downsampled5_reinterpolated5_dsmask_reg.nii.gz', bpref, bpref));
    
%%

% prepare brain --> brainreg registration.
movingVol = double(isonii.img);
movingDims = isonii.hdr.dime.pixdim(2:4);
rMoving = imref3d(size(movingVol), movingDims(2), movingDims(1), movingDims(3));

fixedVol = double(dsregnii.img);
fixedDims = dsregnii.hdr.dime.pixdim(2:4);
rFixed = imref3d(size(fixedVol), fixedDims(2), fixedDims(1), fixedDims(3));

% test to make sure maskVol3 == brainregmaskni.img IT WORKS!
dsmaskregistered = imwarp(dsmasknii.img, rMoving, tform, 'linear', 'OutputView', rFixed);
assert(all(isclose(dsmaskregistered(:), dsmaskregnii.img(:))));

atlVolSize = size(fixedVol);
subVolSize = size(movingVol);

% get the locations in subject volume that correspond to pixel locations in
% the atlas volume 
locVolumeAtlas = getCorrespondingLoc(atlVolSize, tform, rFixed, rMoving, subVolSize);
locVolumeSubject = getCorrespondingLoc(subVolSize, tform.invert, rMoving, rFixed, atlVolSize);

% extract the corresponding patch from the subject volume provided
volRange = arrayfunc(@(x, p) (x: x + p - 1)', location, patchSize);
locPatchAtlas = cellfunc(@(ns) ns(volRange{:}), locVolumeAtlas);


%% compute/load gaussian parameters and estimate unknown pixel values

switch gmmmethod
    case 'compute'
        
        % padding for extracting patches from volume
        locpad = [3 3 3];
        
        if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
        [isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);
        
        % fit gaussian distrubtion (1 cluster)
        gmm = gmdistribution.fit(isopatches,1);
        
        [newPatch, subLocation, newPatchFull] = maximizeRotConditional( gmm.mu, gmm.Sigma, location, locPatchAtlas, movingVol, dsmasknii.img, locVolumeSubject, sigmaMethod, regVal);
        
    case 'load'
        
        wgmmload = load('iso_wgmm_for_katie.mat');
        gmm = wgmmload.isowgmmm3;
        
%         wgmmload = load('wgmm_for_katie.mat');
%         gmm = wgmmload.dswgmmm3;
        
        % extract patch and associated weights
        atlasPatch = fixedVol(volRange{:});
        weightPatch = dsmaskregnii.img(volRange{:});
        weightPatch(weightPatch < 0.001) = 0.001;
        
        % determine the best cluster for the given patch
        logp = gmm.logpost(atlasPatch(:)', weightPatch(:)');
        [~, clust] = max(logp);
        
        [newPatch, subLocation, newPatchFull] = ...
            maximizeRotConditional(gmm.mu(clust,:), gmm.sigma(:,:,clust), location, ...
            locPatchAtlas, movingVol, dsmasknii.img, locVolumeSubject, sigmaMethod, regVal);
        
    case 'groundtruth'
        
        if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
        atlasPatch = niis.iso{subject}.img(volRange{:});
        weightPatch = niis.mask{subject}.img(volRange{:});
        
        % computer ground truth Sigma in atlas space from iso data 
        Sigma = atlasPatch(:)*atlasPatch(:)';
        mu = zeros(1,numel(atlasPatch));
        
%         [newPatch, subLocation, newPatchFull] = ...
%             maximizeRotConditional(mu, Sigma, location, locPatchAtlas, movingVol, dsmasknii.img, ...
%             locVolumeSubject, sigmaMethod, regVal);   
        [newPatch, subLocation, newPatchFull] = ...
            maximizeRotConditional(mu, Sigma, location, locPatchAtlas, double(dsnii.img), dsmasknii.img, ...
            locVolumeSubject, sigmaMethod, regVal);   
        
end

% extract ground truth patch in subject space
volRangeSubj = arrayfunc(@(x, p) (x: x + p - 1)', subLocation, size(newPatch));
realIsoPatch = movingVol(volRangeSubj{:}); 
realIsoPatchwNans = realIsoPatch;
realIsoPatchwNans(isnan(newPatchFull)) = nan;

lininterpPatch = dsnii.img(volRangeSubj{:});
lininterpPatch(isnan(newPatchFull)) = nan;

% visualize
figure(); 
subplot(411); imagesc(subspacetools.reshapeN3Dto2D(realIsoPatchwNans(:)', size(newPatch)), [0 1]); colormap gray; title('Ground Truth');
subplot(412); imagesc(subspacetools.reshapeN3Dto2D(newPatchFull(:)', size(newPatch)), [0 1]); colormap gray; title('Full Patch With Known Pixel Intensities'); 
subplot(413); imagesc(subspacetools.reshapeN3Dto2D(newPatch(:)', size(newPatch)), [0 1]); colormap gray;  title('Patch with Only Estimated Pixel Intensities'); 
subplot(414); imagesc(subspacetools.reshapeN3Dto2D(lininterpPatch(:)', size(newPatch)), [0 1]); colormap gray;  title('Interpolated Patch in Subj Space'); 

%% compare with behaviour of recon
isoPatch = niis.iso{subject}.img(volRange{:});
dsPatch = niis.ds{subject}.img(volRange{:});
weightPatch = niis.mask{subject}.img(volRange{:});
weightPatch(weightPatch < 0.001) = 0.001;
wgmmload = load('iso_wgmm_for_katie.mat');
gmm = wgmmload.isowgmmm3;
rpatches = papago.recon(gmm, dsPatch(:)', weightPatch(:)', 'eig', 95);
rpatches = reshape(rpatches, patchSize);

% visualize
figure(); 
subplot(411); imagesc(subspacetools.reshapeN3Dto2D(isoPatch(:)', patchSize), [0 1]); colormap gray; title('Ground Truth');
subplot(412); imagesc(subspacetools.reshapeN3Dto2D(rpatches(:)', patchSize), [0 1]); colormap gray; title('Full Patch With Known Pixel Intensities'); 
subplot(414); imagesc(subspacetools.reshapeN3Dto2D(dsPatch(:)', patchSize), [0 1]); colormap gray;  title('Interpolated Patch in Subj Space'); 

%%
% 
% brainmoved = imwarp(isonii.img, rMoving, tform, 'linear', 'OutputView', rFixed);
% brainmovedback = imwarp(brainmoved, rFixed, tform.invert, 'linear', 'OutputView', rMoving);
% 
% realIsoPatchRerotated = brainmovedback(volRangeSubj{:});
% % subplot(414); imagesc(subspacetools.reshapeN3Dto2D(realIsoPatchRerotated(:)', size(newPatch)), [0 1]); colormap gray; title('Ground Truth Rotated and Back');
% 
