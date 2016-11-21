function [tgtPatchRange, mins, tgtPatchSize] = srcVol2tgtPatch(srcLoc, srcPatchSize, srcLoc2tgtSpace)
% SRCVOL2TGTPATCH extract the tgt patch range given the src patch location and size.
%
% tgtPatchRange = srcVol2SubjPatch(srcLoc, srcPatchSize, srcLoc2SubjSpace) compute the range of the
% patch in tgt space given the top-left location of the patch in src space. 
%
% srcLoc is a 1-by-nDims vector of the location in src space
% srcPatchSize is a 1-by-nDims vector of the patch size in src apce
% srcLoc2SubjSpace is a 1-by-nDims cell, each entry is of size src
% tgtPatchRange is a 1-by-nDims cell, where each entry is the range of the tgt patch in that
% dimension.
%
% [tgtPatchRange, tgtLoc, tgtPatchSize] = srcVol2tgtPatch(srcLoc, srcPatchSize, srcLoc2SubjSpace)
% also returns the tgt top-left location, and the patch size in tgt space.
% 
%
% TODO: make more general and call from here.

    % get the range of the patch in src space
    srcPatchRange = arrayfunc(@(loc, siz) loc:loc+siz-1, srcLoc, srcPatchSize);
    srcPatchRangeGrid = ndgrid2cell(srcPatchRange{:});
    idx = sub2ind(size(srcLoc2tgtSpace{1}), srcPatchRangeGrid{:});

    % compute the patch limits in tgt space
    mins = cellfun(@(x) floor(min(x(idx(:)))), srcLoc2tgtSpace);
    maxs = cellfun(@(x) ceil(max(x(idx(:)))), srcLoc2tgtSpace);
    
    % compute the final range
    tgtPatchRange = arrayfunc(@(mi, ma) mi:ma, mins, maxs);
    tgtPatchSize = maxs - mins + 1;
    