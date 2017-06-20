function [R, srcMask, tgtMask] = interpmat(srcVolSize, varargin)
%   INTERPMAT return the linear interpolation matrix of a correspondance field

    [R, srcMask, tgtMask] = cor2interpmat(srcVolSize, varargin{:});