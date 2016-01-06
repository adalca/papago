function fwgmm = wgmmfit(X, W, K, varargin)
% WGMMFIT fit a weight - gaussian mixture model
%   fwgmm = wgmmfit(X, W, K)
%   fwgmm = wgmmfit(X, W, K, varargin)
%   
% TODO: comments :)

    fwgmm = wgmm.fit(X, W, K, varargin{:});