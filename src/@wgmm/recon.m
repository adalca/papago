function varargout = recon(wgmm, X, W, method, varargin)
% given weighted gmm, "reconstruct" missing information from data.
%
% methods (varargin):
%   eig (prtile)
%   weig (prtile)
%   pca (prtile)

    error('it seems like we should use svd instead of weig. Unclear, this might've worked from ebginning?');

    inputs = parseInputs(varargin{:});
    varargout = cell(1, nargout);
    
    switch method
        case 'eig'
            varargout{1} = reconeig(wgmm, X, W, inputs.percentile);
           
        case 'weig'
            varargout{1} = reconweig(wgmm, X, W, inputs.percentile);
            
        case 'pca'
            varargout{1} = reconpca(wgmm, X, W, inputs.percentile);
            
        otherwise 
            error('unknown method');
    end
end


function rpatches = reconeig(wgmm, X, W, prctile)
% reconstruct via, for each mixture k:
%   compute principal eigenvectors via eig(sigma(k)), given sigma from WGMM EM
%   project each patch (OP) to this space, 
%   reconstruct patch (RP) via top J scores, where J is chosen to contain prctile of the variance
%   final patch is OP .* w + RP .* (1-w);

    [~, ci] = max(wgmm.logpost(X, W, wgmm.K), [], 2);

    rpatches = zeros(size(X));
    for k = 1:wgmm.K
        cidx = ci == k;
        
        if sum(cidx) == 0
            continue;
        end
        
        % compute eigenvectors and eigenvalues
        [coeff, eigValues] = eig(wgmm.sigma(:, :, k)); % coeff: columns are eigenvectors
        coeff = fliplr(coeff);
        eigValues = fliplr(flipud(eigValues));
        
        % find the number of components that cover at least prctile of variance
        variability = diag(eigValues) ./ sum(diag(eigValues));
        cumvar = cumsum(variability);
        j = find(cumvar > prctile/100, 1, 'first');
        %fprintf('%d components explain %3.1f%% of variance\n', j, cumvar(j)*100);
        
        
        % x difference to cluster center and weights
        mu = wgmm.mu(k, :);
        xc = bsxfun(@minus, X(cidx, :), mu);
        w = W(cidx, :);
        
        % project
        scores = pcax.project(xc', coeff)';
        
        % reconstruct from top j components
        recon = pcax.recon(coeff(:, 1:j), scores(:, 1:j)', mu)';
        rpatches(cidx, :) = w .* X(cidx, :) + (1-w) .* recon;
    end
end

function rpatches = reconweig(wgmm, X, W, prctile)
% reconstruct via, for each mixture k:
%   compute principal eigenvectors via eig(sigma(k)), given sigma from WGMM EM
%   project each patch (OP) to this space, but weight the features of OP and eigenvectors by W.
%   reconstruct patch (RP) via top J scores, where J is chosen to contain prctile of the variance
%   final patch is OP .* w + RP .* (1-w);

    [~, ci] = max(wgmm.logpost(X, W, wgmm.K), [], 2);

    rpatches = zeros(size(X));
    for k = 1:wgmm.K
        cidx = ci == k;
        
        if sum(cidx) == 0
            continue;
        end
        
        % compute eigenvectors and eigenvalues
        [coeff, eigValues] = eig(wgmm.sigma(:, :, k)); % coeff: columns are eigenvectors
        coeff = fliplr(coeff);
        eigValues = fliplr(flipud(eigValues));
        
        % find the number of components that cover at least prctile of variance
        variability = diag(eigValues) ./ sum(diag(eigValues));
        cumvar = cumsum(variability);
        j = find(cumvar > prctile/100, 1, 'first');
        %fprintf('%d components explain %3.1f%% of variance\n', j, cumvar(j)*100);
        
        
        % x difference to cluster center and weights
        mu = wgmm.mu(k, :);
        xc = bsxfun(@minus, X(cidx, :), mu);
        w = W(cidx, :);
        
        % project with weighted components on weighted coefficients
%         scores = zeros(sum(cidx), size(X, 2));
%         for i = 1:sum(cidx)
%             % method 1. 
%             % This is quite slow thouh, since you need to do L' * L and inverse for each patch
%             wc = bsxfun(@times, w(i, :)', coeff); % wc same as (diag * coeff)? Yes, but that is slower
%             wxc = w(i, :)' .* xc(i, :)';
%             scores(i, :) = pcax.project(wxc, wc);
%             
%             % method 2.
%         end
        scores = pcax.project((xc.*W)', coeff)';

        % reconstruct from top j components (note, this time we use coeff, not weighted coeff.
        recon = pcax.recon(coeff(:, 1:j), scores(:, 1:j)', mu)';
        rpatches(cidx, :) = w .* X(cidx, :) + (1-w) .* recon;
    end
end

function rpatches = reconpca(wgmm, X, W, prctile)
% reconstruct via, for each mixture k:
%   compute principal components via PCA of data in this cluster
%   project each patch (OP) to this space, 
%   reconstruct patch (RP) via top J scores, where J is chosen to contain prctile of the variance
%   final patch is OP .* w + RP .* (1-w);

    [~, ci] = max(wgmm.logpost(X, W, wgmm.K), [], 2);

    rpatches = zeros(size(X));
    for k = 1:wgmm.K
        cidx = ci == k;
        
        if sum(cidx) == 0
            continue;
        end
        
        % x difference to cluster center and weights
        mu = wgmm.mu(k, :);
        xc = bsxfun(@minus, X(cidx, :), mu);
        w = W(cidx, :);
        
        % project
        [coeffs, scores, eigValues] = pca(xc(cidx, :));
%         figure(101); imagesc([cov(X(cidx, :)), wgmm.sigma(:,:,k)]);
        
        % find the number of components that cover at least prctile of variance
        variability = eigValues ./ sum(eigValues);
        cumvar = cumsum(variability);
        j = find(cumvar > prctile/100, 1, 'first');

        % reconstruct from top j components
        recon = pcax.recon(coeffs(:, 1:j), scores(:, 1:j)', mu)'; % warning: adding back mu estimated not mu computed
        rpatches(cidx, :) = w .* xc(cidx, :) + (1-w) .* recon;
    end
end


function inputs = parseInputs(varargin)
    
    p = inputParser(); 
    p.addParameter('percentile', 95); % for eigen method
    p.parse(varargin{:});
    inputs = p.Results;
end
