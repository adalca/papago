function [subvolR, srcMask, tgtMask, srcInd, tgtInd, tgtWeight] = ...
    vol2subvolInterpmat(R, srcLoc2tgtSpace, tgtSize, srcSubvolLoc, srcSubvolSize, thr)
% R src --> tgt ===> R is |tgt| x |src|


    if nargin < 6
        % warning('this 0.99 threshold is just emperically set right now...'); 
        thr = 0.99;
    end
    

    % prepare source indeces;
    srcRange = arrayfunc(@(l, s) l:l+s-1, srcSubvolLoc, srcSubvolSize);
    ndSrcRange = ndgrid2cell(srcRange{:});
    srcInd = sub2ind(size(srcLoc2tgtSpace{1}), ndSrcRange{:});

    % prepare target indeces;
    tgtPatchRange = srcVol2tgtPatch(srcSubvolLoc, srcSubvolSize, srcLoc2tgtSpace);
    ndTgtRange = ndgrid2cell(tgtPatchRange{:});
    tgtInd = sub2ind(tgtSize, ndTgtRange{:});
    
    % extract the right region of the R matrix
    subvolR = R(tgtInd, srcInd);
    
    % create the mask of the tgt patch
    % new code
    tgtWeight = reshape(subvolR*ones(numel(srcInd),1), size(tgtInd)) ;
    % old code
    % create the mask of the tgt patch. we do this by forward propogating the target patch and
    % seeing if when you invert the process it returns a value smaller than the original value (in
    % which case it was missing some information from pixels not included)
    % >> tgtWeight = reshape(subvolR*subvolR'*ones(numel(tgtInd),1), size(tgtInd)) ;
    tgtMask = tgtWeight > thr; 
    
    % create the mask of the src patch
    srcMask = ones(size(srcInd)); 
    
    % unnecessary because the source patch will always be smaller than the target patch in this case
    % where we are extracting the patches within the function
    %srcMask = (subvolR'*subvolR*ones(numel(srcInd),1)) > 0.001 ;
    %srcMask = reshape(srcMask, size(srcInd)); 