function [subsampledPatches] = appxDownsizePatches(patches, weights, xlargeSz, largeSz, largeUs, smallUs) 
% [subsampledPatches] = appxDownsizePatches(patches, weights, 13, 9, 5, 1) 

% get blurring filter
warning('whats the best blur size to do? previously was 1/4'); 
sigmaBlur = 1/3 * (largeUs / smallUs);
filtSize = round(sigmaBlur*2.354*2); 
gaussfilter = fspecial3('gaussian',filtSize);

% get the indices corresponding to the inner patch
ConvSz = filtSize + xlargeSz - 1; 
minLoc = floor(ConvSz/2)-floor(largeSz/2);
maxLoc = floor(ConvSz/2)+floor(largeSz/2);
[ii, jj, kk] = ndgrid(minLoc:maxLoc, minLoc:maxLoc, minLoc:maxLoc); 
inds = sub2ind([ConvSz ConvSz ConvSz], ii, jj, kk); 

% create convolution matrix
ConvMtx = convmtxn(gaussfilter, [xlargeSz xlargeSz xlargeSz]);
ConvMtx = ConvMtx(inds,:); 

% blur the patches and reshape
N = size(patches,1); 
blurredPatches = ConvMtx*(patches.*weights)'; 
blurredWeights = ConvMtx*weights'; 
vol4Dstack = reshape(blurredPatches, [largeSz largeSz largeSz, N]);
wt4Dstack = reshape(blurredWeights, [largeSz largeSz largeSz, N]);

% subsample the patches and reshape
subsampledSize = (ceil((ceil(smallUs/largeUs.*largeSz*[1 1 1])/2)+0.5)*2)-1; 
subsampledPatches = volresizeSimple(vol4Dstack./wt4Dstack, [subsampledSize N], 'simplelinear');
%subsampledPatches = volresizeSimple(vol4Dstack./wt4Dstack, [subsampledSize N], 'nearest');
subsampledPatches = reshape(permute(subsampledPatches, [4 1 2 3]), [N prod(subsampledSize)]); 


