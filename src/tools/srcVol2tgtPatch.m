function [tgtPatchRange, mins, tgtPatchSize] = srcVol2tgtPatch(srcLoc, srcPatchSize, srcLoc2tgtSpace)
% SRCVOL2TGTPATCH extract the tgtect patch range given the srcas patch location and size.
%
% tgtPatchRange = srcVol2SubjPatch(srcasLoc, srcPatchSize, srcLoc2SubjSpace) compute the range of the
% patch in tgtect space given the top-left location of the patch in srcas space. 
%
% srcLoc is a 1-by-nDims vector of the location in srcas space
% srcPatchSize is a 1-by-nDims vector of the patch size in srcas apce
% srcLoc2SubjSpace is a 1-by-nDims cell, each entry is of size srcas
% tgtPatchRange is a 1-by-nDims cell, where each entry is the range of the tgtect patch in that
% dimension.
%
% [tgtPatchRange, tgtLoc, tgtPatchSize] = srcVol2tgtPatch(srcasLoc, srcPatchSize, srcLoc2SubjSpace)
% also returns the tgtect top-left location, and the patch size in tgtect space.
% 
%
% TODO: make more general and call from here.

    % get the range of the patch in srcas space
    srcPatchRange = arrayfunc(@(loc, siz) loc:loc+siz-1, srcLoc, srcPatchSize);
    srcPatchRangeGrid = ndgrid2cell(srcPatchRange{:});
    idx = sub2ind(size(srcLoc2tgtSpace{1}), srcPatchRangeGrid{:});

    % compute the patch limits in tgtect space
    mins = cellfun(@(x) floor(min(x(idx(:)))), srcLoc2tgtSpace);
    maxs = cellfun(@(x) ceil(max(x(idx(:)))), srcLoc2tgtSpace);
    
    % compute the final range
    tgtPatchRange = arrayfunc(@(mi, ma) mi:ma, mins, maxs);
    tgtPatchSize = maxs - mins + 1;
    