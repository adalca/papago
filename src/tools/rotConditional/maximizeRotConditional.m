function [newPatch, subLocation, newPatchFull, logP] = ...
    maximizeRotConditional(mu, Sigma, atlLocation, locPatch, subVolume, subWeights, subLocVolume, method, regVal)
% maximize the conditional probablity of a patch's unknown intensities
% given its known intensities in the subject volume by using a specified 
% gaussian model (mu, Sigma).  
% 
% INPUTS
% mu - the mean of the multidimensional Gaussian model
% Sigma - the covariance of the multidimensional Gaussian model
% locPatch - a cell array of patches containing the locations of the
%   corresponding coordinates in the subject volume (subVolume) in double 
%   precision. There is a patch for every dimension (e.g. x, y, z-coordinate)
% subVolume - the subject volume that the patch should be extracted from
% subWeights - a volume of 1's and 0's corresponding to where subVolume 
%   contains known and unknown intensities respectively
%
% OUTPUTS
% newPatch - a patch containing the estimated data in the unknown locations
%   and nans elsewhere
% location - The location corresponding to where the patch is in the subVolume
% newPatchFull - a patch containing the estimated data and known data used

    wmsg = [' the best Using linear interpolation upon rotation. Bilinear probably gives much ', ...
        'better results, though.'];
    sys.warn(wmsg, 'singleWarn', true); 

    % determine the size of the patch
    atlPatchSize = size(locPatch{1}); 

%     % initlize space
%     subPatchSize = zeros(1, nDims);  
%     subLocation = zeros(1, nDims); 
%     modLocPatch = locPatch; 
% 
%     % compute the location and patchsize of the corresponding patch in subject space
%     for i=1:nDims
% 
%         assert(isclean(locPatch{i}), 'Some locations in patch are invalid in subject space'); 
% 
%         minVal = min(floor(locPatch{i}(:))); 
%         maxVal = max(ceil(locPatch{i}(:)));
% 
%         % compute the size of the patch in subject space and the location
%         subPatchSize(i) = maxVal - minVal + 1; 
%         subLocation(i) = minVal; 
% 
%         % modify the locPatch so that the values go from 1 to subPatchSize but
%         % they are double precision
%         modLocPatch{i} = locPatch{i} - minVal + 1; 
%     end
% 
%     % extract the corresponding patch from the subject volume provided
%     volRange = arrayfunc(@(x, p) (x: x + p - 1)', subLocation, subPatchSize);
%     subPatchI = subVolume(volRange{:});
%     subPatchW = subWeights(volRange{:});
% 
%     for i=1:nDims
%         subPatchL{i} = subLocVolume{i}(volRange{:}) - atlLocation(i) + 1; 
%     end

    [modLocPatch, subPatchI, subPatchW, subPatchL, subPatchSize] = ...
        splitSubVols(atlLocation, locPatch, subVolume, subWeights, subLocVolume); 

    
    % calculate the transformation matrix bringing you from subject space to
    % atlas space 
    [R, subMask] = genR(atlPatchSize, subPatchSize, modLocPatch, true(subPatchSize)); 

    switch method
        case 'forward'
            newSigma = pinv(R) * Sigma * pinv(R');
            newSigma = newSigma + eye(size(newSigma)).*regVal;
            %newSigma = nearestSPD(newSigma);
            subMu = pinv(R)*mu(:); 
        case 'inverse'
            [Rinv, subMask] = genR(subPatchSize, atlPatchSize, subPatchL, subMask);
            newSigma = Rinv * Sigma * Rinv';
            subMu = Rinv * mu(:);
        otherwise
            error('not a valid method for sigma estimation');
    end


    subPatchIvalid = subPatchI(subMask(:));
    subPatchWvalid = subPatchW(subMask(:));


    %newSigmaInv = R' * SigmaInv * R;
    %newSigma = pinv(newSigmaInv); 

    % identify which indices we know and which we want to estimate
    unknown = find((subPatchWvalid(:) == 0)); 
    known   = find((subPatchWvalid(:) > 0)); 

    %extract necessary pieces from the newSigma to perform a conditional calculation
    B = newSigma(known, known);
    C = newSigma(known, unknown);

    % compute the data for the missing locations in the new patch
    estimatedData = nan(size(subPatchIvalid)); 
    estimatedData(unknown) = C' * pinv(B) * (subPatchIvalid(known) - subMu(known)); 
    estimatedData(unknown) = estimatedData(unknown) + subMu(unknown); 
    subPatchIvalid(unknown) = estimatedData(unknown); 

    % put the information in a patch of size subPatchSize where there are nan's
    % where we have not calculated any values. 
    newPatch = nan(subPatchSize); 
    newPatch( subMask(:) ) = estimatedData(:); 

    newPatchFull = nan(subPatchSize); 
    newPatchFull( subMask(:) ) = subPatchIvalid(:); 
    
    %compute the marginal probability of being in the current cluster given
    %the known pixels in subject space
    warning('How should we invert B? pinv, inv, \, or add regularization along diagonal?'); 
    logP = -0.5 * (subPatchIvalid(known) - subMu(known))' * pinv(B) *  (subPatchIvalid(known) - subMu(known)); 
    

end

