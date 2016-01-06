function ll = logp(wgmm, varargin)
% compute log likelihood

    logpin = wgmm.logpost(varargin{:}); 
    
    % to void overflow, need to subtract the max in log space. so 
    % sum_n * log (sum_c (pin)); where pin = pi * N()
    % becomes
    % sum_n * log {sum_c (exp(logpin - maxlogpin)) + exp(maxlogpin)}
    %   = sum_n * log {sum_c (exp(logpin - maxlogpin))} + maxlogpin
    maxlogpin = max(logpin, [], 2);
    top = exp(bsxfun(@minus, logpin, maxlogpin));
    sumc = sum(top, 2);
    ll = sum(log(sumc) + maxlogpin, 1);
    
    % check cleanliness
    assert(isclean(ll));
end
