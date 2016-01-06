
% generate the R-transformation matrix using the location information in
% locPatch corresponding to the locations from 1 to subPatchSize
function [R, subMask] = genR(atlPatchSize, subPatchSize, locPatch, subMask)

    nDims = length(locPatch);

    % initlize space
    R = zeros(prod(atlPatchSize), prod(subPatchSize)); 
    location = zeros(nDims, 1); 

    % loop through all locations in the atlas patch and determine which pixels
    % from the subject patch contribute to its intensity
    for i=find(subMask)'

        for j=1:nDims
            pos{j} = [floor(locPatch{j}(i)) ceil(locPatch{j}(i))]; 
            location(j) = locPatch{j}(i); 
        end

        posComb = combvec(pos{:}); 

        if (sum( sum( (posComb > repmat(subPatchSize(:), [1 size(posComb,2)])) | (posComb < 1) ) ) )
            subMask(i) = false;
        else
            q = subvec2ind(subPatchSize, posComb');
            R(i, q) = prod(1 - abs(bsxfun(@minus, posComb, location)));
        end

    end

    % remove regions of the subject patch 
    R(~subMask,:) = []; 

end