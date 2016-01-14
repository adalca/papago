%% initialize
setup

%% parameters
atlPatchSize = [9 9 9];
atlLoc = LOC_LEFT_CORTEX; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;

gmmmethod = 'groundtruth'; % 'compute', 'load', 'groundtruth'
method = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward

subject = 7; %1, 3

%% load volumes
% we'll use some local, small subset, of the buckner data for now.

% TODO: should move this to using the existing md structure, 
% once both katie and adrian have access to the same large data.
bpref = sprintf('buckner%02d', subject);

% subject space
subjisonii = loadNii(sprintf('%s/%s/%s/%s_brain_cropped5.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % iso
dsnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds 
dsmasknii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_dsmask.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));

% atlas space (i.e. registered data)
isoregnii = loadNii(sprintf('%s/%s/%s/%s_brain_iso_2_ds5_us5_size_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds
dsregnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref)); % ds
tf = load(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_reg.mat', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));
dsmaskregnii = loadNii(sprintf('%s/%s/%s/%s_brain_downsampled5_reinterpolated5_dsmask_reg.nii.gz', SYNTHESIS_DATA_PATH, 'buckner/proc', bpref, bpref));
    
%% prepare necessary inputs for reconstruction
% atlLoc, atlPatchSize, subjVol, subjWeightVol, atlLoc2SubjSpace, method
subjVol = double(dsnii.img);
subjWeightVol = logical(dsmasknii.img);

subjDims = dsnii.hdr.dime.pixdim(2:4);
atlDims = dsregnii.hdr.dime.pixdim(2:4);
subjLoc2AtlSpace = tform2cor3d(tf.tform, size(subjVol), subjDims, size(dsregnii.img), atlDims);
atlLoc2SubjSpace = tform2cor3d(tf.tform, size(subjVol), subjDims, size(dsregnii.img), atlDims, 'backward');
extraReconArg = ifelse(strcmp(method, 'inverse'), subjLoc2AtlSpace, regVal);

%% "groundtruth" testing
% atlMu, atlSigma: done below

% get the covariance from the (isotropic) patch we're looking at in atlas space and see if we can
% reconstruction the right patch in subject space.
isoSubjVolInAtl = double(isoregnii.img);

% get mu and sigma in atlas space
atlPatch = cropVolume(isoSubjVolInAtl, atlLoc, atlLoc + atlPatchSize - 1);
atlMu = zeros(atlPatchSize);
atlSigma = atlPatch(:) * atlPatch(:)' + eye(prod(atlPatchSize)) * regVal;

% reconstruct
[reconPatch, logp, subjPatchLoc] = paffine.recon(atlMu, atlSigma, atlLoc, atlPatchSize, subjVol, subjWeightVol, atlLoc2SubjSpace, method, extraReconArg);

% look at "true" patch
initSubjPatch = cropVolume(dsmasknii.img, subjPatchLoc, subjPatchLoc + size(reconPatch) - 1);
initSubjPatch(isnan(reconPatch)) = nan;
correctSubjPatch = cropVolume(subjisonii.img, subjPatchLoc, subjPatchLoc + size(reconPatch) - 1);
correctSubjPatch(isnan(reconPatch)) = nan;

view3Dopt(correctSubjPatch, reconPatch, initSubjPatch)

% quick visualization
subplot(1, 3, 1); imagesc(patchview.reshapeto2D(correctSubjPatch)); 
subplot(1, 3, 2); imagesc(patchview.reshapeto2D(reconPatch));
subplot(1, 3, 3); imagesc(patchview.reshapeto2D(initSubjPatch));
colormap gray;





%% compute/load gaussian parameters and estimate unknown pixel values

switch gmmmethod
    case 'compute' % compute the covariance from a patch column
        
%         % padding for extracting patches from volume
%         locpad = [3 3 3];
%         
%         if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
%         [isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, atlPatchSize, atlLoc, locpad);
%         
%         % fit gaussian distrubtion (1 cluster)
%         gmm = gmdistribution.fit(isopatches,1);
%         
%         [newPatch, subLocation, newPatchFull] = maximizeRotConditional( gmm.mu, gmm.Sigma, atlLoc, locPatchAtlas, movingVol, dsmasknii.img, subjLoc2AtlSpace, method, regVal);
%         
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
            maximizeRotConditional(gmm.mu(clust,:), gmm.sigma(:,:,clust), atlLoc, ...
            locPatchAtlas, movingVol, dsmasknii.img, subjLoc2AtlSpace, method, regVal);
        
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
            maximizeRotConditional(mu, Sigma, atlLoc, locPatchAtlas, double(dsnii.img), dsmasknii.img, ...
            subjLoc2AtlSpace, method, regVal);   
        
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
rpatches = reshape(rpatches, atlPatchSize);

% visualize
figure(); 
subplot(411); imagesc(subspacetools.reshapeN3Dto2D(isoPatch(:)', atlPatchSize), [0 1]); colormap gray; title('Ground Truth');
subplot(412); imagesc(subspacetools.reshapeN3Dto2D(rpatches(:)', atlPatchSize), [0 1]); colormap gray; title('Full Patch With Known Pixel Intensities'); 
subplot(414); imagesc(subspacetools.reshapeN3Dto2D(dsPatch(:)', atlPatchSize), [0 1]); colormap gray;  title('Interpolated Patch in Subj Space'); 

%%
% 
% brainmoved = imwarp(isonii.img, rMoving, tform, 'linear', 'OutputView', rFixed);
% brainmovedback = imwarp(brainmoved, rFixed, tform.invert, 'linear', 'OutputView', rMoving);
% 
% realIsoPatchRerotated = brainmovedback(volRangeSubj{:});
% % subplot(414); imagesc(subspacetools.reshapeN3Dto2D(realIsoPatchRerotated(:)', size(newPatch)), [0 1]); colormap gray; title('Ground Truth Rotated and Back');
% 
