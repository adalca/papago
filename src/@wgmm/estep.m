function [ll, expect, wg] = estep(wg, data)
% e-step 
%
% since this is a mixture model, most model variants will include an update of cluster membership
% probabilities, usually indicated by a variable gamma. They are usually some form of 
%   gamma <-- p(clust) * p(data|clust) / sum_clust(p(clust)*p(data|clust).
% where the data includes some weighting aspect.
%
% depending on the model, the e-step might include other expectation updates. 

    % get (unnormalized) posterior
    % this is log( p(clust) * p(data|clust) )
    logpin = wg.logpost(data); % N x k
    
    % to void overflow, need to subtract the max in log space. so 
    % gamma = pin / sum(pin) where pin = pi * N()
    % becomes
    % gamma = exp(logpin - maxlogpin + maxlogpin) / sum(exp(logpin - maxlogpin + maxlogpin))
    %       = exp(logpin - maxlogpin) * exp(maxlogpin) / sum(exp(logpin - maxlogpin) * exp(maxlogpin))
    %       = exp(logpin - maxlogpin) / sum(exp(logpin - maxlogpin))
    maxlogpin = max(logpin, [], 2);
    top = exp(bsxfun(@minus, logpin, maxlogpin));
    ll = sum(log(sum(top, 2)) + maxlogpin, 1);
    assert(isclean(ll), 'log likelihood is not clean');
    
    % compute expected cluster assignment
    expect.gammank = bsxfun(@rdivide, top, sum(top, 2));
    assert(isclean(expect.gammank), 'cluster assignments are not clean');
    
    % output the number of points in each cluster
    if wg.opts.verbose >= 2
        mi = argmax(expect.gammank, [], 2);
        clustAsnHist = hist(mi, 1:data.K);
        nDigits = numel(num2str(max(clustAsnHist)));
        fprintf('Clust Mem: ');
        fprintf(sprintf('%%%dd', nDigits+1), clustAsnHist); fprintf('\n');
    end
end
