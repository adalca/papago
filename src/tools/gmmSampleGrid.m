function [samples, csamples] = gmmSampleGrid(gmm, patchSize, gridSize, sigma)
% sample from a gmm object, and arrange samples in a grid image for ease of display
%   samples = gmmSampleGrid(gmm, patchSize, gridSize)
%   samples = gmmSampleGrid(gmm, patchSize, gridSize, sigma)
%   [samples, csamples] = gmmSampleGrid(...)
%
% gmm: a gmdistribution obj
% patchSize: a 1 x nDims vector of the size of the patches. 
% gridSize: a 1 x 2 vector of the size of the grid in which to arrange samples
% sigma (optional) is the same size as gmm.Sigma 
%
% samples is a nClust cell array, each element a grid of samples
% csamples is a nClust cell array, each element a matrix of vectorized samples, nSamples x nDims

    % parse inputs
    nClust = size(gmm.mu, 1); % is there a cleaner way to do this?
    assert(size(gmm.mu, 2) == prod(patchSize));
    if nargin <= 3
        sigma = gmm.Sigma;
    end
    
    % sample
    samples = cell(nClust, 1);
    csamples = cell(nClust, 1);
    for c = 1:nClust
        csamples{c} = mvnrnd(gmm.mu(c, :), sigma(:,:,c), prod(gridSize)); 
        
        % split into a cell array of patches
        csamplescell = dimsplit(1, csamples{c});
        csamplescell = cellfunc(@(x) reshape(x, patchSize), csamplescell);
        csamplescell = reshape(csamplescell, gridSize);
        
        % combine into a large array
        samples{c} = cell2mat(csamplescell);
    end
    