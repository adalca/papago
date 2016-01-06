function S = posterior2assignment(post, forcefullassignment)
% use the posterior points in a gmm to decide on the point assignment to clusters. 
% forcefullassignment forces each cluster to have at least one element, which we achieve
% heuristically by going over each emtpy cluster (after initial assignment) and grabbing the most
% likely points for that cluster among points from "large" clusters.
%
% not checked too thoroughly :)
        
        [~, S] = max(post, [], 2);
        K = size(post, 2);

        if forcefullassignment
            % for those clusters that have zero votes, give them some reasonable points
            h = hist(S, 1:K); 
            z = find(h == 0);
            allbigclusters = find(h > (numel(z)+1));
            ptsinbigclusters = find(any(bsxfun(@eq, S, allbigclusters), 2));
            for i = 1:numel(z), 
                k = z(i); 
                [~, pti] = max(post(ptsinbigclusters, k)); 
                S(ptsinbigclusters(pti)) = k; 
                ptsinbigclusters(pti) = [];
            end
            h = hist(S, 1:K); 
            assert(all(h > 0));
        end
        