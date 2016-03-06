function [cropArray] = mdBoundingBoxRange(md, inputModality, outputModality, varargin)
    % bounding boxes

    mods = {inputModality};
    if nargin > 2 && ~isempty(outputModality)
        mods = {inputModality, outputModality};
    end
    
    [cropArray] = md.applyfun(@boundingBoxRangeNii, mods, varargin{:});
end
