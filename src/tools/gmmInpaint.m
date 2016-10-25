function Xi = gmmInpaint(Xi, gmd, ci, clusts)
% gmmInpaint
%   
% Xi = gmmInpaint(Xi, gmd) inpaint data given MATLAB's gmdist object
%
% TODO: move to papago or wgmm?

    % parse inputs and prepare variables.
    nClust = numel(gmd.ComponentProportion);
    if nargin <= 2 || isempty(ci)
        % get cluster assignment
        [p, ci] = max(gmd.posterior(Xi), [], 2); 
        assert(isclean(p));
    end
    
    if nargin <= 3
        clusts = 1:nClust;
    end
    


    % inpaint
    for c = clusts
        cidx = find(ci == c);
        
        Xc = bsxfun(@minus, Xi(cidx, :), gmd.mu(c, :));
        
        % loop because pcax.inpaint only works on one nan pattern at a time.
        for k = 1:numel(cidx) 
            Xi(cidx(k), :) = pcax.inpaint(Xc(k, :)', gmd.Sigma(:,:,c))';
        end
        Xi(cidx, :) = bsxfun(@plus, Xi(cidx, :), gmd.mu(c, :));
    end
    