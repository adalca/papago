function [subjPatchRange, mins, subjPatchSize] = atl2SubjPatch(atlasLoc, atlPatchSize, atlLoc2SubjSpace)
% ATL2SUBJPATCH extract the subject patch range given the atlas patch location and size.
%
% subjPatchRange = atl2SubjPatch(atlasLoc, atlPatchSize, atlLoc2SubjSpace) compute the range of the
% patch in subject space given the top-left location of the patch in atlas space. 
%
% atlasLoc is a 1-by-nDims vector of the location in atlas space
% atlPatchSize is a 1-by-nDims vector of the patch size in atlas apce
% atlLoc2SubjSpace is a 1-by-nDims cell, each entry is of size atlas
% subjPatchRange is a 1-by-nDims cell, where each entry is the range of the subject patch in that
% dimension.
%
% [subjPatchRange, subjLoc, subjPatchSize] = atl2SubjPatch(atlasLoc, atlPatchSize, atlLoc2SubjSpace)
% also returns the subject top-left location, and the patch size in subject space.
% 
%
% TODO: make more general and call from here.

    % get the range of the patch in atlas space
    atlasPatchRange = arrayfunc(@(loc, siz) loc:loc+siz-1, atlasLoc, atlPatchSize);
    atlasPatchRangeGrid = ndgrid2cell(atlasPatchRange{:});
    idx = sub2ind(size(atlLoc2SubjSpace{1}), atlasPatchRangeGrid{:});

    % compute the patch limits in subject space
    mins = floor(min(cellfun(@(x) x(idx(:)), atlLoc2SubjSpace)));
    maxs = ceil(max(cellfun(@(x) x(idx(:)), atlLoc2SubjSpace)));
    
    % compute the final range
    subjPatchRange = arrayfunc(@(mi, ma) mi:ma, mins, maxs);
    subjPatchSize = ma-mi+1;
    