function datasetLocHist(atlvol, volumes, varargin)
% DATASETLOCHIST Show a histogram at a location accross a (rigidly registered) volume dataset
%   datasetLocHist(atlvol, volumes) for 2D data
%   datasetLocHist(atlvol, volumes, slice) for 3D data
%
% atvol is the 3D atlas
% volumes is a 4D volume: the 3D subject volumes stacked in the 4th dimension. 

    patches = reshape(volumes, [numel(atlvol), 1, size(volumes, 4)]);
    [~, grididx] = patchlib.vol2lib(atlvol, ones(1, ndims(atlvol)));
    patchview.votehist(atlvol, patches, grididx, ones(1, ndims(atlvol)), varargin{:});