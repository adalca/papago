function varargout = recon(wg, X, W, method, varargin)
% given weighted gmm, "reconstruct" missing information from data.
%
% methods (varargin):
%   eig (prtile)
%   weig (prtile)
%   pca (prtile)

    %error('it seems like we should use svd instead of weig. Unclear, this might''ve worked from ebginning?');

    inputs = parseInputs(varargin{:});
    varargout = cell(1, nargout);
    
    switch method
        case 'eig'
            varargout{1} = reconeig(wg, X, W, inputs.percentile);
           
        case 'weig'
            varargout{1} = reconweig(wg, X, W, inputs.percentile);
            
        case 'pca'
            varargout{1} = reconpca(wg, X, W, inputs.percentile);
            
        case 'latentMissing'
    
            xRecon = X;
            obsIdx = W == 1;
            name = wg.opts.model.name;
            wg.opts.model.name = 'latentMissing';
            [~, expect] = wg.estep(struct('Y', X, 'W', W));
            maxk = argmax(expect.gammank, [], 2);
            wg.opts.model.name = name;
            
            X(~obsIdx) = nan;

            for i = 1:size(X, 1)
                lobsIdx = obsIdx(i, :);
                k = maxk(i);

                oosigmak = wg.params.sigma(lobsIdx, lobsIdx, k);
                mosigmak = wg.params.sigma(~lobsIdx, lobsIdx, k);
                b = (X(i, lobsIdx) - wg.params.mu(k, lobsIdx))';
                t1 = mosigmak * (oosigmak \ b); 
                xRecon(i, ~lobsIdx) = wg.params.mu(k, ~lobsIdx) + t1';
                
                % test using inpaintWithGaussConditional. This is the same, but much slower due to
                % asserts and overhead.
                % xReconCheck = inpaintWithGaussConditional(X(i, :)', wgmm.params.mu(k,:)', ...
                %    wgmm.params.sigma(:, :, k), true);
                % assert(all(isclose(xRecon(i, :), xReconCheck')));
            end
            varargout{1} = xRecon;
            varargout{2} = expect;
            
        case 'latentSubspace'
            K = size(wg.expect.gammank, 2);
            Y = X;
            
            [dHigh, dLow, ~] = size(wg.params.W);
            
            % compute Xhat
            muk = zeros(K, size(Y,2));
            denom = zeros(K, size(Y,2));
            S = zeros(dLow, dLow, size(Y,1), K);
            Xhat = zeros(size(Y,1), dLow, K);
            for i = 1:size(Y, 1)
                obsIdx = W(i, :) == 1;
                yobs = Y(i, obsIdx);

                for k = 1:K
                    sigmasq_t = wg.params.sigmasq(k);
                    w = wg.params.W(obsIdx, :, k);
                    L = w ./ sqrt(sigmasq_t);
                    
                    % ppca() uses Sherman-Morrison, probably because it's safer?
                    % S_ki = inv(L * L' + eye(dLow));
                    ltl = L'*L;
                    S_ki = eye(dLow) - ltl / (eye(dLow) + ltl);
                    X_ki = S_ki/sigmasq_t * (w' * (yobs - wg.params.mu(k, obsIdx))');
                    
                    % update mu_k
                    muterm = yobs' - w * X_ki;
                    muk(k, obsIdx) = muk(k, obsIdx) + wg.expect.gammank(i, k) .* muterm';
                    
                    % denominator for mu_k
                    denom(k, obsIdx) = denom(k, obsIdx) + wg.expect.gammank(i, k);
                    
                    % update large matrices
                    S(:,:,i,k) = S_ki;
                    Xhat(i,:,k) = X_ki';
                end
            end
            
            xRecon = X*0;
            xReconChk = X*0;
            mi = argmax(wg.expect.gammank, [], 2);
            
            for k = 1:K
                w = wg.params.W(:,:,k);
                
                wwt = w'*w;
                cRecon = w / (wwt) * (wwt + wg.params.sigmasq(k) * eye(dLow)) * Xhat(mi==k, :, k)';
                xRecon(mi==k, :) = bsxfun(@plus, cRecon', wg.params.mu(k, :));
                
                xReconChk(mi==k, :) = bsxfun(@plus,  Xhat(mi==k, :, k) * w', wg.params.mu(k, :));
            end
            varargout{1} = xRecon;
            varargout{2} = xReconChk;
            
        case 'latentMissingR'

            xRecon = wg.opts.model.data.yorig;
            [~, expect] = wg.estep(struct('Y', X, 'W', W));

            maxk = argmax(expect.gammank, [], 2);
            
            R = wg.opts.model.data.R;
            ydsmasks = wg.opts.model.data.ydsmasks;
            yorigs = wg.opts.model.data.yorig;
            
            sigmac = dimsplit(3, wg.params.sigma);
            
            for i = 1:size(X, 1)
                obsIdx = ydsmasks{i};
                yobs = yorigs{i}(obsIdx);
                r = R.data(R.idx{i}, :); 
                k = maxk(i);
                
                sigmaz = sigmac{k};
                
                % rotate sigma. see paffine.atl2SubjGauss
                sigmai = r * sigmaz * r';
                mui = wg.params.mu(k, :) * r';
                
                % extract sigmas
                oosigmak = sigmai(obsIdx, obsIdx);
                mosigmak = sigmai(~obsIdx, obsIdx);
                
                y_ik = yorigs{i};
                y_ik(~obsIdx) = mui(~obsIdx) + ...
                    (mosigmak * (oosigmak \ (yobs - mui(obsIdx))'))';
                xRecon{i} = y_ik;
            end
            varargout{1} = xRecon;

            
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
