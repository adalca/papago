function [nanpatches, obsIdx] = destroyFeatures(patches, nDestroy, method)
% destory features (voxels) in patches. Destroyed features will appear as NaNs
%
% patches is nSamples x prod(patchSize)
% nDestory is scalar specifying how many voxels to destroy
% method is 'rand-consistent' or 'rand-inconsistent'. 
%   rand-inconsistent means destroy voxels at random, but the same voxels for each patch
%   rand-inconsistent means destroy voxels at random independent for each patch
%   TODO: destroy in the plane sense, the way it happens in real data.
%
% nanpatches are the destroyed patches
%
% patch subspace project

    % prepare inputs and useful variables.
    narginchk(3, 3);
    assert(isclean(patches));
    [nSamples, nElems] = size(patches);
    invalidvox = [true(1, nDestroy), false(1, nElems - nDestroy)];
    
    % initalize nanpatches
    nanpatches = patches;

    % destroy features/voxels
    switch method
        case 'rand-consistent'
            % ranomize which voxels are invalid
            invalid = invalidvox(randperm(nElems));
            nanpatches(:, invalid) = nan;
            
        case 'rand-inconsistent'
            % ranomize which voxels are invalid, but differently for different patches
            invalid = arrayfunc(@(x) invalidvox(randperm(nElems)), 1:nSamples);
            invalid = cat(1, invalid{:});
            nanpatches(invalid) = nan;
            
        otherwise
            error('destroyFeatures: unknown method');
    end
    
    obsIdx = ~invalid;
    