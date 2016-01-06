function locVolume = getCorrespondingLoc(atlVolSize, tform, rAtlas, rSubject, subVolSize)
% get the locations in subject volume that correspond to pixel locations in
% the atlas volume given the specified affine transformation tform and
% return the locations in a cell array where each cell corresponds to a
% dimension
%
% Inputs:
%
% altVolSize - the size of the atlas volume. We will find a corresponding
%   location in the subject volume for each point in the atlas volume
% tform - the transformation that relates subject space and atlas space. It
%   was found using imregtform. 
% rAtlas - information about the pixel sizes in atlas space, obtained using
%   imref3d
% rSubject - information about the pixel sizes in subject space, obtained
%   using imref3d
% subVolSize - the size of the subject volume. If the corresponding point
%   in atlas space is outside of this size, we put a NaN in that location
%
% Output:
% 
% locVolume - a cell array containing information about the location in
%   the subject volume that corresponds to the associated pixel in atlas
%   space. For instance locVolume{1}(10,5,3) gives us the index of the 1st
%   dimension of the corresponding point in the subject volume that point
%   (10,5,3) in atlas space corresponds to. 
%   e.g. subVolume(loc{1}(10,5,3), loc{2}(10,5,3), loc{3}(10,5,3)) ~=  atlVolume(10,5,3)

    % get the number of dimensions, number of voxels, and size of the atlas volume
    nDims       = length(atlVolSize); 
    nVoxels     = prod(atlVolSize); 

    % generate vectors to use in ndgrid2cell
    for i=1:nDims
        ndgridVec{i} = 1:atlVolSize(i); 
    end

    % get matricies that give the coordinate of each pixel in cell structure
    atlPosCell = ndgrid2cell(ndgridVec); 

    % vectorize each cell and add it as a row in the matrix atlPos
    atlPosCellVec = cellfunc(@(x) x(:), atlPosCell);

    % compute the corresponding position in the subject volume 
    %rAtlasPos = rAtlas.intrinsicToWorld(atlPosCellVec{:});
    %rSubjectPos = transformPointsInverse(tform, rAtlasPos);
    %subPos = rSubject.worldToIntrinsic(rSubjectPos);

    subPos = zeros(3,nVoxels); 
    [x, y, z] = rAtlas.intrinsicToWorld(atlPosCellVec{2}, atlPosCellVec{1}, atlPosCellVec{3});
    [xx, yy, zz] = transformPointsInverse(tform, x, y, z);
    [subPos(2,:), subPos(1,:), subPos(3,:)] = rSubject.worldToIntrinsic(xx, yy, zz);

    %atlPos = cat(2, atlPosCellVec{:}); 
    %subPos = tform' \ [atlPos'; ones(1, nVoxels)]; 

    % make sure you are within the valid range of the subject volume
    subPos(subPos < 1) = nan; 
    subPos(subPos > repmat(subVolSize(:), [1 nVoxels])) = nan; 

    % reshape the volume and place into cell array
    locVolume = dimsplit(1, subPos); 
    locVolume = cellfun(@(x) reshape(x, atlVolSize), locVolume, 'UniformOutput', false); 

