function spangrid = showImageSpace(imgdata, x, nSpace, shape, varargin)
% compute and show representative images in the space spanned by the given coordinates x
%
% imgdata is N x D
% x is N x 1 or N x 2
% nSpace is N x 1 or N x 2
% shape is 1 x 2

    % parse inputs
    [x, nSpace, opts] = parseInputs(imgdata, x, nSpace, shape, varargin{:});
    
    % span the space.
    xi = span(x, nSpace);
    
    % compute the normalized weights
    wt = weights(xi, x, opts.lambda);
    
    % compute the weighted sums.
    spanvec = squeeze(sum(bsxfun(@times, wt, permute(imgdata, [3, 1, 2])), 2));
    
    % plot the images
    plotImages(spanvec, xi, shape, opts.plotArgs{:}, 'nSpace', nSpace);
   
    % compute a grid of images 
    spangrid = datavec2grid(spanvec, shape, nSpace);
end

function wt = weights(xi, x, lambda)
% compute weights via wt = exp(-lambda * (xi - x) * S^{-1} * (xi - x)'

    S = cov(x) + eps;
    dst = pdist2(xi, x, 'mahalanobis', S);
    wt = exp(-lambda * dst .* 2);
    wt = bsxfun(@times, wt, 1./sum(wt, 2));
end

function xi = span(x, nSpace)
% span the space sampled by x. x is N x 1 or N x 2

    xi = cellfunc(@(x, n) linspace(min(x), max(x), n), dimsplit(2, x), dimsplit(2, nSpace(:)'));
    xi = ndgrid2cell(xi{:});
    xi = cellfunc(@(x) x(:), xi); % vectorize
    xi = cat(2, xi{:}); % cell to mat
    
    assert(size(xi, 1) == prod(nSpace));
end

function spanvalimg = datavec2grid(spanval, shape, nSpace)
% take data from vector form (N x D) to a grid (nSpace .* shape)
% e.g. a 100 x 25 could go to a 10x10 grid of 5x5 images.

    spanvalc = dimsplit(1, spanval);
    spanvalc = cellfunc(@(x) reshape(x, shape), spanvalc);
    spanvalc = reshape(spanvalc, nSpace);
    spanvalimg = cell2mat(spanvalc);
end

function [x, nSpace, opts] = parseInputs(data, x, nSpace, shape, varargin)
% parse inputs.

    p = inputParser();
    p.addRequired('data', @ismatrix);
    p.addRequired('x', @ismatrix);
    p.addRequired('nSpace', @(n) size(x, 2) == numel(n));
    p.addRequired('shape', @(n) all(size(n) == [1, 2]));
    p.addParameter('lambda', 5, @isscalar);
    p.addParameter('plotArgs', {'border', 1, 'borderColor', [1, 0, 0]}, @iscell);
    p.parse(data, x, nSpace, shape, varargin{:});
    opts = p.Results;
    
    % transform a 2
    if size(x, 2) == 1
        x = [x, zeros(size(x, 1), 1)];
        nSpace = [nSpace, 1];
    end
end
