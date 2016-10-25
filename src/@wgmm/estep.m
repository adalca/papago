function [ll, gammank, varargout] = estep(wgmm, X, W)
% e-step 
%
% since this is a mixture model, most model variants will include an update of cluster membership
% probabilities, usually indicated by a variable gamma. They are usually some form of 
%   gamma <-- p(clust) * p(data|clust) / sum_clust(p(clust)*p(data|clust).
% where the data includes some weighting aspect.
%
% depending on the model, the e-step might include other expectation updates. 

    varargout = cell(0);
    if strcmp(wgmm.logpUpdateMethod, 'model5')
        varargout = cell(1);
    end

    % get (unnormalized) posterior
    % this is log( p(clust) * p(data|clust) )
    [logpin, varargout{:}] = wgmm.logpost(X, W); % N x k
    
    % to void overflow, need to subtract the max in log space. so 
    % gamma = pin / sum(pin) where pin = pi * N()
    % becomes
    % gamma = exp(logpin - maxlogpin + maxlogpin) / sum(exp(logpin - maxlogpin + maxlogpin))
    %       = exp(logpin - maxlogpin) * exp(maxlogpin) / sum(exp(logpin - maxlogpin) * exp(maxlogpin))
    %       = exp(logpin - maxlogpin) / sum(exp(logpin - maxlogpin))
    maxlogpin = max(logpin, [], 2);
    top = exp(bsxfun(@minus, logpin, maxlogpin));
    ll = sum(log(sum(top, 2)) + maxlogpin, 1);
    
    % compute expected cluster assignment
    wgmm.expect.gammank = bsxfun(@rdivide, top, sum(top, 2));
    assert(isclean(wgmm.expect.gammank));
    
    % TODO: only do this if have any sum(gammank, 2) == 0
    % Should change! should destory clusters that don''t have any support.
    % warning('estep: adding eps to gamma. TODO: Fix this!');
    % gammank = gammank + eps; 
    
    % for certain models, we have other useful E-step quantities
    if strcmp(wgmm.logpUpdateMethod, 'model5')
        wgmm.expect.Xk = varargout{1};
    end
    
    % check cleanliness of variable
    assert(isclean(ll));
end