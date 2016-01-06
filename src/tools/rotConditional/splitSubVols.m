function [modLocPatch, subPatchI, subPatchW, subPatchL, subPatchSize] = ...
    splitSubVols(atlLocation, locPatch, subVolume, subWeights, subLocVolume)
% maximize the conditional probablity of a patch's unknown intensities
% given its known intensities in the subject volume by using a specified
% gaussian model (mu, Sigma).
%
% INPUTS
% locPatch - a cell array of patches containing the locations of the
%   corresponding coordinates in the subject volume (subVolume) in double
%   precision. There is a patch for every dimension (e.g. x, y, z-coordinate)
% subVolume - the subject volume that the patch should be extracted from
% subWeights - a volume of 1's and 0's corresponding to where subVolume
%   contains known and unknown intensities respectively
%
% OUTPUTS

% determine the number of dimensions for the patch and the size of the patch
nDims = length(locPatch);

% initlize space
subPatchSize = zeros(1, nDims);
subLocation = zeros(1, nDims);
modLocPatch = locPatch;

% compute the location and patchsize of the corresponding patch in subject space
for i=1:nDims
    
    assert(isclean(locPatch{i}), 'Some locations in patch are invalid in subject space');
    
    minVal = min(floor(locPatch{i}(:)));
    maxVal = max(ceil(locPatch{i}(:)));
    
    % compute the size of the patch in subject space and the location
    subPatchSize(i) = maxVal - minVal + 1;
    subLocation(i) = minVal;
    
    % modify the locPatch so that the values go from 1 to subPatchSize but
    % they are double precision
    modLocPatch{i} = locPatch{i} - minVal + 1;
end

% extract the corresponding patch from the subject volume provided
volRange = arrayfunc(@(x, p) (x: x + p - 1)', subLocation, subPatchSize);
subPatchI = subVolume(volRange{:});
subPatchW = subWeights(volRange{:});

for i=1:nDims
    subPatchL{i} = subLocVolume{i}(volRange{:}) - atlLocation(i) + 1;
end