function [R, srcMask, tgtMask] = cor2interpmat(srcVolSize, varargin)
% COR2INTERPMAT return the linear interpolation matrix of a correspondance field
% 
% [R, srcMask, tgtMask] = cor2interpmat(srcVolSize, tgtLoc2SrcSpace) return the (sparse)
% interpolation matrix R of a correspondance field, i.e. that takes the source and warps is to the
% target. In other words, tgtVoxels = R * srcVoxes, up to linear interpolation
% 
% The function requires the correspondance field tgtLoc2SrcSpace, in double precision, of each
% target point to the original source volume e.g. such as returned by tform2cor3d(tform, ...,
% 'backward'). tgtLoc2SrcSpace is a 1-by-nDims cell, where each entry is a volume of the same size
% as the target. the target size is then not required as an input. tgtLoc2SrcSpace{d}(i) is then the
% dth coordinate of the ith input. srcVolSize is a 1-by-nDims vector for the size of the source
% volume.
%
% The returned R matrix is tgtVoxels-by-srcVoxels. srcMask is a srcVolSize boolean volume, meaning
% which voxels *contribute* to the target. tgtMask is a tgtVolSize boolean volume, and refers to
% which voxels were *fully covered* by source voxels. Note that some rows and columns of R can be
% all zeros -- to strip R of these rows or columns simply run 
%   >> R(~tgtMask, ~srcMask) = [];
%
% [R, srcMask, tgtMask] = cor2interpmat(srcVolSize, tgtVolSize, tgtLoc2SrcSpaceVec) allows
% tgtLoc2SrcSpaceVec to be specified in matrix format, where tgtLoc2SrcSpaceVec is
% prod(tgtVolSize)-by-nDims.
%
% Note: design decision: if *any* of the source point contributing to target point i is outside the
% scope of the source, target point i becomes nan. This is a strict decision, another options might
% to still estimate target point i from the available source points, even if not all the source
% points needed are available. This is less acurate, but you get more points estimated, especially
% around the edges. TODO: perhaps we should make this an option?
% 
% TOTEST: test interpolation with this and compare with volwarp (linear).
%
% part of papago project. We mostly use this for patch correspondances between two spaces (e.g.
% subject space and atlas space). 

    % parse inputs. tgtLoc2SrcSpace will be a cell of size 1-by-nDims
    [srcVolSize, tgtVolSize, tgtLoc2SrcSpace] = parseInputs(srcVolSize, varargin{:});
    nDims = numel(srcVolSize);
    
    % initialize transform-interpolation matrix
    srcMask = false(srcVolSize);
    tgtMask = false(tgtVolSize);
    
    % get floor and ceil of tgtLoc2SrcSpace, temporary but useful matrices below
    floorTgtLoc2SrcSpace = cellfunc(@(x) floor(x), tgtLoc2SrcSpace);
    ceilTgtLoc2SrcSpace = cellfunc(@(x) ceil(x), tgtLoc2SrcSpace);
    
    % prepare R sparse matrix constructs
    Rivec = zeros(2^nDims * prod(tgtVolSize), 1);
    Rjvec = zeros(2^nDims * prod(tgtVolSize), 1);
    Rvec = zeros(2^nDims * prod(tgtVolSize), 1);
    
    % go through each voxel of the target, and for each voxel get the 2^nDims source points that
    % will contribute linearly to this target point.
    for tgti = 1:prod(tgtVolSize)
        
        % get the floor and ceil of tgtLoc2SrcSpace(i) in each dimension. In other words, these
        % points determine the hypercube bordering the double-precision position in the source
        % tgtLoc2SrcSpace(i) corresponding to ith target voxel. gridLim is 1-by-nDims, each entry is
        % 1-by-2.
        gridLim = cellfunc(@(f, c) [f(tgti), c(tgti)], floorTgtLoc2SrcSpace, ceilTgtLoc2SrcSpace);
        
        % make sure we don't have an actual integer value, and if we do, only use that point once.
        % Note: we shouldn't just use unique() as that is very slow over so many calls. Check first
        % if needed.
        isEqual = cellfun(@(x) x(1) == x(2), gridLim);
        if any(isEqual)
            gridLim = cellfunc(@(x) unique(x), gridLim);
        end
        
        % get all the grid points
        gridPosCell = cellfunc(@(x) x(:), ndgrid2cell(gridLim{:}));
        gridPos = cat(2, gridPosCell{:});
        
        % check if all the grid points are actually within the source volume
        % if all grid poitns are 
        if all(gridPos >= 1 & bsxfun(@le, gridPos, srcVolSize))
            srcMask(gridPosCell{:}) = true; % this means that pixel *contributed*
            tgtMask(tgti) = true;
            
            % compute the relative source cube position compared to the actual (double-precision)
            % source position tgtLoc2SrcSpace(i). 
            tgtLoc2SrcSpacei = cellfun(@(x) x(tgti), tgtLoc2SrcSpace);
            relativeGridPos = bsxfun(@minus, gridPos, tgtLoc2SrcSpacei(:)'); % (2^nDims)-by-nDims
            
            % write R: for each source point on the grid, take the product of (1 - dst-to-pt)
            srci = sub2ind(srcVolSize, gridPosCell{:}); 
            
            % update sparse R constructs
            idx = 2^nDims * (tgti-1) + (1:2^nDims);
            idx = idx(1:numel(srci));
            Rivec(idx) = tgti;
            Rjvec(idx) = srci;
            Rvec(idx) = prod(1 - abs(relativeGridPos), 2); % linear interpolation
        end
    end
    
    % construct the sparse matrix
    Rvec(Rivec==0) = [];
    Rjvec(Rivec==0) = [];
    Rivec(Rivec==0) = [];
    R = sparse(Rivec, Rjvec, Rvec, prod(tgtVolSize), prod(srcVolSize));
end

function [srcVolSize, tgtVolSize, tgtLoc2SrcSpace] = parseInputs(srcVolSize, varargin)
    
    % checkn number of inputs.
    narginchk(2, 3);
    
    % interpmat(srcVolSize, tgtLoc2SrcSpace)
    if nargin == 2
        tgtLoc2SrcSpace = varargin{1};
        tgtVolSize = size(tgtLoc2SrcSpace{1});
        
    % interpmat(srcVolSize, tgtVolSize, tgtLoc2SrcSpaceVec)
    else
        tgtVolSize = varargin{1};
        tgtLoc2SrcSpaceVec = varargin{2};
        
        % check inputs
        assert(size(tgtLoc2SrcSpaceVec, 1) == prod(tgtVolSize), ...
            'Number of target voxels do not match')
        assert(size(tgtLoc2SrcSpaceVec, 2) == numel(srcVolSize), ...
            'Dimensions of tgtLoc2SrcSpaceVec does not match tgtVolSize');
        
        % get tgtLoc2SrcSpace
        tgtLoc2SrcSpaceCell = dimsplit(2, tgtLoc2SrcSpaceVec);
        tgtLoc2SrcSpace = cellfunc(@(x) reshape(x, varargin{1}), tgtLoc2SrcSpaceCell);
    end
    
    % check that volumes have the same # of dimensions
    assert(numel(tgtVolSize) == numel(srcVolSize), 'Number of volume dimensions do not agree');
end
