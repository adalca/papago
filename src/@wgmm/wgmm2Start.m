function S = wgmm2Start(gmm, starttype, varargin)
% WGMM2Start creates a Start structure from a wgmm or gmdistribution object
%   S = wgmm2Start(gmm, starttype, ...) creates a Start structure from a wgmm or gmdistribution
%   object. starttype can be 'parameters', in which case S is a struct with 'ComponentProportion',
%   'Sigma', and 'mu', or 'assignment', in which case wgmm2Start also receives the posterior data.
%
% See Also: wgmm, posterior2assignment 

switch starttype
    case 'assignment'
        post = varargin{1};
        S = wgmm.posterior2assignment(post, true);

    case 'parameters'
        if isa(gmm, 'gmdistribution')
            S.ComponentProportion = gmm.ComponentProportion;
            S.Sigma = gmm.Sigma;
            S.mu = gmm.mu;
            
        else
            assert(isa(gmm, 'wgmm'));
            S.ComponentProportion = gmm.pi;
            S.Sigma = gmm.sigma;
            S.mu = gmm.mu;
        end
        
    otherwise
        error('wgmm2Start: Unknown start type');
end