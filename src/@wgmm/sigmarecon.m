function sigmar = sigmarecon(sigma, wtw, method)
% SIGMARECON reconstruct weakly voted covariance matrix sigma
%
% sigmar = sigmarecon(sigma, wtw, method) compute sigma reconstruction (due to areas of low weight)
% each entry sigma_ij is "reconstructed" somehow using entries s_ik and s_jk, etc. Sigma is a square
% covariance matrix, wtw is a square weight matrix of the same size as sigma indicating the amount
% of weight (votes) used to build sigma (wtw = W' * W; where W is nPts x nFeat). method can be:
%   'weighted-mean-1'
%   'greedy-1'
%
% part of the wgmm project.
%
% contact: {adalca,klbouman}@csail.mit.edu
    
    switch method
        case 'wmean1'
            % use a weighted sum over the correlations going through all x_k
            % specifically, s(i, j) = wmean(w, q)
            % where w = w_ik * w_jk, and q = s_ik * s_jk / s_kk
            %
            % sigmareg = sigma * 0;
            % for i = 1:size(sigma, 1)
            %     for j = 1:size(sigma, 2)
            %         q = sigma(i, :) .* sigma(j, :) ./ diag(sigma)';
            %         w = wtw(i, :) .* wtw(j, :);
            % 
            %         % weighted mean can be one of several:
            %         %  s2(i, j) = wmean(w, q);
            %         sigmareg(i, j) = wmean(w .* q, q);
            %     end
            % end

            
            % trying to speed regularization up.
            d = diag(sigma)';
            s1 = bsxfun(@times, sigma, 1./d);
            allsisj = permute(bsxfun(@times, s1, permute(sigma, [3, 2, 1])), [1, 3, 2]);
            allw = permute(bsxfun(@times, wtw, permute(wtw, [3, 2, 1])), [1, 3, 2]);
            sigmar = sum(allw .* allsisj .* allsisj, 3) ./ sum(allw .* abs(allsisj), 3);
        
        case 'greedy1'
            % the greedy-1 method takes a maximum weighted correlated pixel x_k, and uses the
            % correlation of x_k with x_j as an indication for the correlation of x_i with x_j.
            % The weighted correlation used is sigma_ik * w_ik * w_jk / sqrt(sigma_kk).
            % For final sigmareg is then s_ik * s_jk / s_kk
            %
            % % slow (double-loop) implementation:
            % sigmar = sigma * 0;
            % for i = 1:size(sigma, 1)
            %     for j = 1:size(sigma, 2) 
            %         teststat = sigma(i, :) .* wtw(i, :) .* wtw(j, :) ./ sqrt(diag(sigma)');
            %         [~, mi] = max(teststat);
            %         sigmar(i, j) = sigma(i, mi) .* sigma(j, mi) ./ sigma(mi, mi);
            %     end
            % end
            
            
            % single-loop implementation. 
            % This is a good medium between no-loop (memory intense) and double-loop (slow)
            d = diag(sigma)';
            sigmar = sigma * 0;
            z = bsxfun(@times, sigma .* wtw, 1 ./ sqrt(diag(sigma)'));
            for i = 1:size(sigma, 1)
                teststat = bsxfun(@times, z(i, :), wtw);
                [~, mi] = max(teststat, [], 2);
                
                x = ndgrid2cell(i, 1:size(sigma, 1));
                ii = sub2ind(size(sigma), x{1}(:), mi(:));
                ij = sub2ind(size(sigma), x{2}(:), mi(:));
                
                sigmar(i, :) = sigma(ii) .* sigma(ij) ./ d(mi)';
            end
            
            % % fast implementation. but this is very memory-intense.
            % d = diag(sigma)';
            % dsq = sqrt(d);
            % allw = permute(bsxfun(@times, wtw, permute(wtw, [3, 2, 1])), [1, 3, 2]);
            % m = bsxfun(@rdivide, bsxfun(@times, allw, permute(sigma, [1, 3, 2])), permute(dsq, [1, 3, 2]));
            % [~, mi] = max(m, [], 3);
            % 
            % x = ndgrid2cell(1:size(sigma, 1), 1:size(sigma, 1));
            % ii = sub2ind(size(sigma), x{1}(:), mi(:));
            % ij = sub2ind(size(sigma), x{2}(:), mi(:));
            % 
            % sigmar = reshape(sigma(ii) .* sigma(ij) ./ d(mi(:))', size(sigma));

        case 'none'
            sigmar = zeros(size(sigma));
            
        otherwise
            error('wgmm.sigmareg: Unknown method');
    end
end
