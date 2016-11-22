function [ll, expect, wgmm] = estep(wgmm, data)
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
    logpin = wgmm.logpost(data); % N x k
    
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
    % TODO: >> expect.gammank = expect.gammank + eps; 
    %   need to consider whether to add eps to gammank. 
    %   or only do this if have any sum(gammank, 2) == 0
    %   or should destroy clusters that don't have any support
    [expect, wgmm] = reclusterHeuristic(wgmm, expect);
    assert(isclean(expect.gammank), 'cluster assignments are not clean');
end

function [expect, wg] = reclusterHeuristic(wg, expect)

    if isfield(wg.opts.model, 'reclusterThreshold')
        thr = wg.opts.model.reclusterThreshold;
    else
        thr = 100;
    end

    mi = argmax(expect.gammank, [], 2);
    h = hist(mi, 1:size(expect.gammank, 2));
    
    if any(h < thr)
        fprintf('%6d', h); fprintf('\n');
        fprintf(2, 'E-Step: %d clusters have less than %d datapoints\n%s\n', ...
            sum(h < thr), thr, 'using heuristics to re-define clusters.');
        
        largestCluster = argmax(h);
        % option 1: force a few random ! datapoints to be part of the new cluster, and set the
        % parameters to the overall mean/covar? This implementation isn't finished, since here we
        % were just grabbing datapoints from the largest cluster and not re-setting the
        % parameters.
        %             idx = find(mi == largestCluster);
        %             ridx = idx(randperm(numel(idx), 100));
        %             mi(ridx) = i;
        %             o = zeros(1, numel(h));
        %             o(i) = 1;
        %             expect.gammank(ridx, :) = repmat(o, [100, 1]);
        
        % option 2: split the largest cluster into two.
        f = find(h < thr);
        for fi = 1:numel(f)
            i = f(fi);
            idx = find(mi == largestCluster);
            nIdx = numel(idx);
            nSelIdx = floor(nIdx/2);
            assert(nSelIdx > thr, 'Don''t have enough data to resample');
            selidx = idx(randsample(nIdx, nSelIdx));
            
            % force new points to be in the new cluster
            o = zeros(1, numel(h));
            o(i) = 1;
            expect.gammank(selidx, :) = repmat(o, [nSelIdx, 1]);
            h(i) = h(i) + nSelIdx;
            
            % reset parameters. This might be complicated due to the different modes.
            wg.params.mu(i,:) = wg.params.mu(largestCluster, :);
            wg.params.sigma(:,:,i) = wg.params.sigma(:, :,largestCluster);
            wg.params.pi(i) = wg.params.pi(largestCluster);
            if isfield(wg.params, 'W')
                wg.params.W(:,:,i) = wg.params.W(:, :,largestCluster);
            end
            if isfield(wg.params, 'sigmasq')
                wg.params.sigmasq(i) = wg.params.sigmasq(largestCluster);
            end
            
            mi = argmax(expect.gammank, [], 2);
            h = hist(mi, 1:size(expect.gammank, 2));
            largestCluster = argmax(h);
        end
        
        h = hist(mi, 1:size(expect.gammank, 2));
        fprintf('%6d', h); fprintf('\n');
    end
    
end
