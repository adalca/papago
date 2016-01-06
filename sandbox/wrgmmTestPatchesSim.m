%% setup
subspacesetup; 
patchSize = [9 9 9];
location = LOC_LEFT_CORTEX; 
locpad = [3 3 3]; 

%% extract patches
if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
[dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
[maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
[isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);

%% learn a GMM and sample patches from that GMM 

N = 1000;                   % number of patches to sample from GMM
D = size(isopatches, 2);    % dimensionality of the patch
K = 1;                      % number of clusters in GMM to learn from data

% learn and sample from a gmm from the nLearn random patches.
nlearn = min(size(isopatches, 1), max(10000, D)+1);
r = randperm(size(isopatches, 1));
Xlearn = isopatches(r(1:nlearn), :);
[X, mu, sigma, piDist, cidx, sigmainv] = subspacetools.simgmm(N, D, K, 'frompatchesgmm', Xlearn);

% make the mean 0
X = bsxfun(@minus, X, mu(:)'); 
mu = zeros(size(mu)); 

% put the sampled patches along with weights into weighted gmm class
W = ones(size(X)); 
wgmmsampled = wgmm(X, W, K, mu, sigma, piDist, sigmainv); 

% visualize the means of the learned GMM 
figure(); imagesc(subspacetools.reshapeN3Dto2D(wgmmsampled.mu, patchSize)); colormap gray; title('means'); axis equal off;

% visualize some new samples vs real patches
subspacetools.comparePatchGMMData(Xlearn, ones(nlearn, 1), X, ones(N, 1), patchSize, 30);
subplot(1,2,1); title('real data'); subplot(1,2,2); title('resampled data');

%% 

% generate a bunch of 2D rotation angles in degrees
rotAngle = 90*(2*rand(N,1) - 1); %70*ones(N,1); 
%rotAngle = 70 * ones(N, 1); 

% calculate the size of the rotated patches
largePatchSize = [ceil(2*sqrt(patchSize(1:2).^2/2)) 1];

% initilize space to save the patches, weights and Rotation matricies
largePatches = zeros(N, prod(largePatchSize));
smallPatches = zeros(N, prod(patchSize));
weights = zeros(N, prod(patchSize));
Rmtx = zeros(prod(patchSize), prod(largePatchSize), N);

% generate subject space weight mask
subjectMask = zeros(largePatchSize);

offset = 3; 
spacingOn = 1:offset:largePatchSize(1);
spacingOff = 1:largePatchSize(1);
spacingOff(spacingOn) = []; 

subjectMask(spacingOn, :) = 1;

% identify which indices we know and which we want to estimate
unknown = find(subjectMask(:) == 0); 
known = find(subjectMask(:) > 0 ); 
idxOrder = [unknown; known];

% compute the centers of the large and small patches for rotation
center1 = (patchSize(1:2)+1)./2;
center2 = (largePatchSize(1:2)+1)./2;

% iterate through all of the sampled patches
for p=1:N
    
    % extract a small patch and rotate it using bilinear interpolation
    smallPatch = reshape(X(p,:), patchSize);
    largePatch = imrotate(smallPatch, -rotAngle(p), 'bilinear');
    largePatches(p,:) = largePatch(:);
    
    % compute the rotation matrix that brings largePatch back to ~smallPatch
    count = 1;
    for j=1:patchSize(2)
        for i=1:patchSize(1)
            
            % get pixel location in small patch and subtract small patch center
            ii = i - center1(1);
            jj = j - center1(2);
            
            % compute 2D rotation matrix
            theta = -rotAngle(p)*pi/180;
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            
            % rotate pixel location and add back the center of the large patch
            outPos = R * [ii; jj];
            outPos = outPos + center2';
            
            % determine the 4 pixels to take information from
            outPos1 = [floor(outPos(1)); floor(outPos(2))];
            outPos2 = [floor(outPos(1)); ceil(outPos(2))];
            outPos3 = [ceil(outPos(1)); ceil(outPos(2))];
            outPos4 = [ceil(outPos(1)); floor(outPos(2))];
            
            % put bilinear weights in rotation matrix 
            Rmtx(count, sub2ind(largePatchSize, outPos1(1), outPos1(2) ), p) = prod(1 - abs(outPos1 - outPos));
            Rmtx(count, sub2ind(largePatchSize, outPos2(1), outPos2(2) ), p) = prod(1 - abs(outPos2 - outPos));
            Rmtx(count, sub2ind(largePatchSize, outPos3(1), outPos3(2) ), p) = prod(1 - abs(outPos3 - outPos));
            Rmtx(count, sub2ind(largePatchSize, outPos4(1), outPos4(2) ), p) = prod(1 - abs(outPos4 - outPos));
            
            count = count + 1;
            
        end
    end
    
    % using the determined rotation matrix re-compute the small patch and weights
    largePatch = largePatch.*subjectMask; 
    
    % linearly interpolate missing values
    [Xpos, Ypos] = meshgrid(1:largePatchSize(1), 1:largePatchSize(2)); 
    largePatch2  = largePatch; 
    largePatch2(spacingOff,:) = interp2(Xpos(spacingOn, :),Ypos(spacingOn, :),largePatch(spacingOn, :),Xpos(spacingOff, :),Ypos(spacingOff, :));
   
    smallPatches(p,:) = Rmtx(:,:,p)*largePatch2(:);
    weights(p,:) = Rmtx(:,:,p)*subjectMask(:); 
      
end

%%

warning('small patches are wrong because I didnt delete information in the large patches'); 
 
weights(weights<0.000001) = 0.000001;

% since the smallPatches will be a little different than the original small
% patches then we recalculate the GMM covariance
wgmmObj = wgmmfit(smallPatches, weights, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 10);

rpatches = papago.recon(wgmmObj, smallPatches, weights, 'pca', 95); 
papago.visRecon(smallPatches, weights, {rpatches}, patchSize, 'titles', {'recon'},  'nShow', 3);


warning('take this out - just for debugging'); 
weightsold = weights; 
weights = ones(size(weights)); 
wgmmObj = wgmmfit(smallPatches, weights, K, 'debug', false, 'sigmareg', 0.00001, 'replicates', 10);



%NOTE MEAN MUST BE 0 FOR THIS
figure;
for p=1:N
    
    %sigmaNewInv = Rmtx(:,:,p)' * diag(weights(p,:)) * wgmmObj.sigmainv * diag(weights(p,:)) * Rmtx(:,:,p);
    %sigmaNew = pinv(sigmaNewInv); 

    
    % we want to compute the new sigma, to do this we have to invert what I
    % had previously written down
    sigmaNew = pinv( Rmtx(:,:,p) ) * diag(1./weights(p,:)) * wgmmObj.sigma * diag(1./weights(p,:)) * pinv( Rmtx(:,:,p)' );
    sigmaNew = sigmaNew + eye(size(sigmaNew)).*1e-5;
    
%     J = sigmaNewInv([unknown; known], [unknown; known]); 
%     
%     J11 = sigmaNewInv(unknown,unknown);
%     J12 = sigmaNewInv(unknown,known); 
%     h1 = -J12*largePatches(p,known)'; 
%     
%     newPatch1 = nan(largePatchSize);
%     newPatch1(unknown) = pinv(J11)*(h1 - J12*largePatches(p,known)'); 
%     newPatch1(known) = largePatches(p,known);

    
    %extract necessary pieces from the sigmaNew to perform a condiional calculation
    B = sigmaNew(known, known);
    C = sigmaNew(known, unknown);
    
    % compute the data for the missing locations in the new patch
    newPatch2 = nan(largePatchSize);
    newPatch2(unknown) = C' * pinv(B) * largePatches(p,known)'; 
    newPatch2(known) = largePatches(p,known);
    
    % rotate back the rotated patch
    rotPatch = reshape(rpatches(p,:), patchSize);
    largeRotPatch = imrotate(rotPatch, -rotAngle(p), 'bilinear');

    
    largePatch = reshape(largePatches(p,:), largePatchSize); 
    imagesc([largeRotPatch largePatch newPatch2], [-.5 .5]); axis equal; axis off; colormap gray; 
    %subplot(122); imagesc([largePatch([2 4 6 8 10 12], :) newPatch2([2 4 6 8 10 12],:) ]); axis equal; axis off; colormap gray;
    %imwrite([largeRotPatch largePatch newPatch2] + 0.5, sprintf('/Users/klbouman/Downloads/result_%d.png', p)); 
    pause(2);
end
