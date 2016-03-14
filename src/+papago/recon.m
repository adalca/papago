function varargout = recon(wgmm, X, W, method, varargin)
% given weighted gmm, "reconstruct" missing information in patch considering the papago setting.
%
% methods (varargin):
%   cond () % X and W can be cell?

    inputs = parseInputs(varargin{:});
    varargout = cell(1, nargout);
    
    switch method
        case 'cond' % conditional reconstruction
            [varargout{1:2}] = cond(wgmm, X, W, inputs);
            
        otherwise 
            error('unknown method');
    end
end


function [rpatches, rlocs] = cond(wgmm, X, W, inputs)
% TODO: maximizeRotConditional should not be the one that also returns the location?
    
    logp = wgmm.logpost(X, W);
    [~, clust] = max(logp);

    % Point estimate
    [~, rlocs, rpatches] = ...
        maximizeRotConditional(wgmm.mu(clust,:), wgmm.sigma(:,:,clust), ...
        inputs.location, inputs.locPatchAtlas, ...
        inputs.subjdsvol, inputs.subjmask, ...
        inputs.locVolumeSubject, inputs.reconSigmaMethod, inputs.reconRegVal);
    
%     % check in-plane 
%     cr = arrayfunc(@(r, p) r:(r + p - 1), reconSubs{i}{t}, size(reconPatches{i}{t}));
%     maps = ~isnan(reconPatches{i}{t}) & subjmask(cr{:});
%     q = subjdsvol(cr{:});
%     assert(all(isclose(reconPatches{i}{t}(maps), q(maps))));
end



function inputs = parseInputs(varargin)
    
    p = inputParser(); 
    
    % for rot estimates
    p.addParameter('location', [], @isvector);
    p.addParameter('locPatchAtlas', []);
    p.addParameter('locVolumeSubject', []);
    p.addParameter('subjdsvol', []);
    p.addParameter('subjmask', []);
    p.addParameter('reconSigmaMethod', 'inverse'); % inverse of forward
    p.addParameter('reconRegVal', 1e-4, @isscalar);
    p.parse(varargin{:});
    
    inputs = p.Results;

end