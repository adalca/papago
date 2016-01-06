function [gammank, ll] = estep(wgmm, X, W)
% e-step (posteriors)

    % get posterior
    logpin = wgmm.logpost(X, W); % N x k
    
    % to void overflow, need to subtract the max in log space. so 
    % gamma = pin / sum(pin) where pin = pi * N()
    % becomes
    % gamma = exp(logpin - maxlogpin + maxlogpin) / sum(exp(logpin - maxlogpin + maxlogpin))
    %       = exp(logpin - maxlogpin) * exp(maxlogpin) / sum(exp(logpin - maxlogpin) * exp(maxlogpin))
    %       = exp(logpin - maxlogpin) / sum(exp(logpin - maxlogpin))
    maxlogpin = max(logpin, [], 2);
    top = exp(bsxfun(@minus, logpin, maxlogpin));
    gammank = bsxfun(@rdivide, top, sum(top, 2));
    ll = sum(log(sum(top, 2)) + maxlogpin, 1);
    
    % TODO: only do this if have any sum(gammank, 2) == 0
    % Should change! should destory clusters that don''t have any support.
    % warning('estep: adding eps to gamma. TODO: Fix this!');
    % gammank = gammank + eps; 
    
    % check cleanliness of variable
    assert(isclean(gammank));
    assert(isclean(ll));
end