function wg = recluster(wg)
% recluster wg if necessary.
%
% TODO: >> expect.gammank = expect.gammank + eps; 
%   need to consider whether to add eps to gammank. 
%   or only do this if have any sum(gammank, 2) == 0
%   or should destroy clusters that don't have any support

    % extract important features
    thr = wg.opts.reclusterThreshold;
    method = wg.opts.reclusterMethod;

    % check for small clusters
    mi = argmax(wg.expect.gammank, [], 2);
    clustAsnHist = hist(mi, 1:size(wg.expect.gammank, 2));
    
    if any(clustAsnHist < thr)
        fprintf('%6d', clustAsnHist); fprintf('\n');
        fprintf(2, 'E-Step: %d cluster(s) have less than %d datapoints\n%s\n', ...
            sum(clustAsnHist < thr), thr, 'using heuristics to re-define clusters.');
        
        switch method
            case 'splitLargest'

                % we'll need to know the largest cluster
                largestCluster = argmax(clustAsnHist);
                
                % (repeatedly, if necessary) split the largest cluster into two
                lowClusters = find(clustAsnHist < thr);
                for fi = 1:numel(lowClusters)
                    lowClustIdx = lowClusters(fi);
                    
                    % get points that (still) belong to the largest cluster
                    idx = find(mi == largestCluster);
                    nIdx = numel(idx);
                    nSelIdx = floor(nIdx/2); % select half of the points
                    assert(nSelIdx > thr, 'Don''t have enough data to resample');
                    selidx = idx(randsample(nIdx, nSelIdx));
                    
                    % force new points to be in the new cluster
                    selGammank = zeros(1, numel(clustAsnHist));
                    selGammank(lowClustIdx) = 1;
                    wg.expect.gammank(selidx, :) = repmat(selGammank, [nSelIdx, 1]);
                    
                    % reset parameters
                    wg.params.mu(lowClustIdx, :) = wg.params.mu(largestCluster, :);
                    if isfield(wg.params, 'sigma')
                        wg.params.sigma(:, :, lowClustIdx) = wg.params.sigma(:, :, largestCluster);
                    end
                    wg.params.pi(lowClustIdx) = wg.params.pi(largestCluster)/2;
                    wg.params.pi(largestCluster) = wg.params.pi(largestCluster)/2;
                    % also update low-dim parameters
                    if isfield(wg.opts.model.dopca, 'W')
                        wg.params.W(:, :, lowClustIdx) = wg.params.W(:, :, largestCluster);
                        wg.params.sigmasq(lowClustIdx) = wg.params.sigmasq(largestCluster);
                    end
                    
                    % recompute histogram of cluster assignment and largest cluster
                    mi = argmax(wg.expect.gammank, [], 2);
                    clustAsnHist = hist(mi, 1:size(wg.expect.gammank, 2));
                    largestCluster = argmax(clustAsnHist);
                end
            case 'addEps'
                wg.expect.gammank = wg.expect.gammank + eps; % does not add up to 1
            
            otherwise
                error('Unknown Reclustering Method');
        end
        
        fprintf('%6d', clustAsnHist); fprintf('\n');
    end    
end
