function vols = dctnBasis(inputSize, nBases)
% e.g. inputSize = 5,5,5, nBases = 3,3,3 (up to 3 in e/a dimension
% following http://www-rohan.sdsu.edu/doc/matlab/toolbox/images/transfo6.html

    narginchk(1,2)
    if nargin == 1
        nBases = inputSize;
    end

    ndims = numel(inputSize);
    assert(ndims == numel(nBases));

    % get all the points of execution
    x = arrayfunc(@(x) 1:x, nBases);
    ndvec = ndgrid2vec(x{:});
%     cvec = ndgrid2cell(x{:});
%     idx = sub2ind(nBases, cvec{:});
    
    % go through each volume i need to create
    vols = cell(nBases);
    for i = 1:size(ndvec,1)
        
        % get ndgrid of nd voxels
        x = arrayfunc(@(x) 0:(x-1), inputSize);
        imvec = ndgrid2vec(x{:});
        
        cosArg = bsxfun(@times, pi * (2 * imvec + 1) / 2, (ndvec(i,:)-1)./inputSize);
        im = prod(cos(cosArg), 2) * prod(alpha(ndvec(i,:), inputSize), 2);
        vols{i} = reshape(im, inputSize);
    end
end

function x = alpha(p, m)
    assert(all(p > 0));
    x = (p == 1) .* sqrt(2./m) + (p > 1) .* sqrt(2./m);
end