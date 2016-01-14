function srcLoc2TgtSpace = tform2cor3d(tform, srcVolSize, srcVoxDims, tgtVolSize, tgtVoxDims, dirn)
% TFORM2COR transform a 3D tform to a coordinate correspondance matrix
%   
% srcLoc2TgtSpace = tform2cor(tform, srcVolSize, srcDims, tgtVolSize, tgtDims) transform a 3D affine
% tform (affine3d object) between source and target volumes, to a coordinate correspondance matrix.
% tform can be obtained using something like imregtform(source, ..., target, ...). srcVolSize and
% tgtVolSize are the sizes of the source and target volumes, respectively. srcVoxDims and tgtVoxDims
% are the voxel dimensions of the source and target (e.g. 1x1x7mm ==> [1, 1, 7]). srcLoc2TgtSpace is
% a cell array of 1-by-nDims, where each entry is the same size as the source volume (srcVolSize).
% Each voxel in these entries represents the (double-precision) coordinate of source voxel in target
% space according to the tform. 
%
% cor = tform2cor(tform, srcVolSize, srcDims, tgtVolSize, tgtDims, dirn) where dirn is a boolean. if
% dirn is 'forward', cor = srcLoc2TgtSpace. if dirn is 'backward', cor is the backwards
% correspondance matrix, i.e. tgtLoc2SrcSpace, whose cell entries are of the size of the target;
% Note, in both cases tform is the same transform: the one that takes the source to the target.
%
% Part of papago project
%
% Example:
% assuming we have a tform obtained using 
% >> tform = imtform(subj, subjRef, atlas, atlasRef, ...)
% we can do:
% >> subjLoc2AtlSpace = tform2cor3d(tform, size(subj), subjVoxDims, size(atlas), atlasVoxDims);
% >> atlLoc2SubjSpace = tform2cor3d(tform.invert, size(atlas), atlasVoxDims, size(subj), subjVoxDims);

    % check inputs and set invert default to false.
    narginchk(5, 6);
    if nargin == 5, dirn = 'forward'; end
    
    % if requiring
    if strcmp(dirn, 'forward')
        % prepare a vectorize nd-grid of the source volume size
        srcPos = cellfunc(@(x) x(:), size2ndgrid(srcVolSize));

        % create coordinate reference objects (imref3d) for subject and target
        % these are necessary structures to take us from 
        srcRef = imref3d(srcVolSize, srcVoxDims(2), srcVoxDims(1), srcVoxDims(3));
        tgtRef = imref3d(tgtVolSize, tgtVoxDims(2), tgtVoxDims(1), tgtVoxDims(3));

        % go from source intrinsic coordinates to target intrinsic coordinates:
        [wx, wy, wz] = srcRef.intrinsicToWorld(srcPos{2}, srcPos{1}, srcPos{3});
        [xx, yy, zz] = transformPointsForward(tform, wx, wy, wz);
        srcLoc2TgtSpace = zeros(prod(srcVolSize), 3);
        [srcLoc2TgtSpace(:, 2), srcLoc2TgtSpace(:, 1), srcLoc2TgtSpace(:, 3)]  = ...
            tgtRef.worldToIntrinsic(xx, yy, zz);

        % take out any rows with coordinates outside of the volume
        srcLoc2TgtSpace(srcLoc2TgtSpace < 1 | bsxfun(@gt, srcLoc2TgtSpace, tgtVolSize)) = nan; 
        srcLoc2TgtSpace(any(isnan(srcLoc2TgtSpace), 2), :) = nan; 

        % reshape the volume and place into cell array
        srcLoc2TgtSpace = cellfunc(@(x) reshape(x, srcVolSize), dimsplit(2, srcLoc2TgtSpace));
        
    else
        assert(strcmp(dirn, 'backward'));
        srcLoc2TgtSpace = tform2cor3d(tform.invert, tgtVolSize, tgtVoxDims, srcVolSize, srcVoxDims);
    end

    assert(all(cellfun(@(x) all(x(:) >= 1 | isnan(x(:))), srcLoc2TgtSpace)), ...
        'output has coordinates beyond limits');