function [gmms, patches] = trainhuh(volumefile, idx)
% train (weighted) gaussian mixture models given volumes. 
%
%
% operation: given a column of subvolumes, 
%   1. load the necessary patches
%   2. do training (learn gmms)
%   3. return group of GMMs
%
% Warning: location has to be relative to this stack of subvolumes
    
    nDims = numel(atlVolSize);

    % setup main gmm grid
    assert(all(gridSpacing > 0));
    subvolOverlap = subvolSize - gridSpacing;
    [subgridlocs, subGridLength] = papago.prepgrid(atlVolSize, subvolSize, subvolOverlap);


    volumes



    
    
    % build the necessary gmm or wgmm
    gmms = cell(1, numel(gmmMethods));
    for t = 1:numel(gmmMethods)
        
        % train gmm
        gmms{t} = papago.train(gmmMethods{t}, patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
    end


    warning('trying subtraction of means');
    
    
    K, gmmargs, wgmmargs, reconSubj, varargin
    
    
    
    
    
    gmm = learngmm(patches, weights, K, opt, wgmmargs, varargin)