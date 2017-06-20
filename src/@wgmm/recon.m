function varargout = recon(wg, data, method, varargin)
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
            varargout{1} = reconeig(wg, Y, W, inputs.percentile);
           
        case 'weig'
            varargout{1} = reconweig(wg, Y, W, inputs.percentile);
            
        case 'pca'
            varargout{1} = reconpca(wg, Y, W, inputs.percentile);
            
        case 'latentMissing'
            if wg.opts.verbose > 0
                warning('LM reconstruction currenty uses most likely cluster, not wsum')
            end
            
            % extract the data
            Y = data.Y;
            if isfield(data, 'W')
                obsMask = data.W;
                assert(islogical(obsMask) | all(obsMask(:) == 0 | obsMask(:) == 1));
                obsMask = obsMask == 1; % in case it's not logical
            else
                obsMask = ~isnan(Y);
            end
            N = size(Y, 1);
            
            % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than indexing
            % into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times per iteration, it
            % ends up mattering for our iterations.
            sigmaPrevCell = dimsplit(3, wg.params.sigma);
            
            % get maximum cluster assignment via e-step. Make sure the Estep is done in LM
            if strcmp(wg.opts.model.name, 'latentMissingR')
                modelname = wg.opts.model.name;
                wg.opts.model.name = 'latentMissing';
            end
            [~, expect] = wg.estep(data);
            maxk = argmax(expect.gammank, [], 2);
            if strcmp(wg.opts.model.name, 'latentMissingR')
                wg.opts.model.name = modelname;
            end
            
            % prepare output variables
            yRecon = Y;
            for i = 1:N
                % extract the observed y indices
                obsIdx = obsMask(i, :);
                yobs = Y(i, obsIdx);
                
                % most likely cluster
                k = maxk(i);

                % extract sigma and mu
                sigmaPrev = sigmaPrevCell{k};
                oosigmak = sigmaPrev(obsIdx, obsIdx);
                mosigmak = sigmaPrev(~obsIdx, obsIdx);
                mu = wg.params.mu(k, :);
                
                % reconstruct y.
                yRecon(i, ~obsIdx) = mu(~obsIdx) + (mosigmak * (oosigmak \ (yobs - mu(obsIdx))'))';
                
                1;
            end
            varargout{1} = yRecon;
            varargout{2} = expect.gammank;

        case 'latentSubspaceR'
            
            misThr = 0.97;
            
            % data
            R = data.R;                     % R rotations from atlas to subject space.
            Y = data.Y;                     % data in subject space
            ydsmasks = data.ydsmasks;       % down-sample mask cell -- (should be "planes" in subject space)
            ydsmasksFullVoxels = data.ydsmasksFullVoxels;
            N = numel(Y);
            assert(data.K == size(wg.params.mu,1));
            [~, dLow, ~] = size(wg.params.W);
            
            % get maximum cluster assignment via e-step. Make sure the Estep is done in LMR
            name = wg.opts.model.name;
            wg.opts.model.name = 'latentMissingR';
            [lll, expect] = wg.estep(data);
            maxk = argmax(expect.gammank, [], 2);
            wg.opts.model.name = name;
            
            % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows -- so we
            % transpose the R data matrix and index in columns instead of rows.
            RdataTrans = R.data';
            
            rcondwarnings = 0;
            
            % compute Xhat
            yRecon = cell(1, N);
            yReconChk = cell(1, N);
            for i = 1:N
                k = maxk(i);
                
                % get the data
                obsIdx = ydsmasksFullVoxels{i};
                ySubj = Y{i};
                ySubjObs = ySubj(obsIdx);
                
                % compute statistics
                r = RdataTrans(:, R.idx{i})';
                v = wg.params.sigmasq(k);
                w = r * wg.params.W(:, :, k);
                wobs = w(obsIdx, :);
                L = wobs ./ sqrt(v);
                muSubj = r * wg.params.mu(k, :)';

                % ppca() uses Sherman-Morrison, probably because it's safer?
                % S_ki = inv(L * L' + eye(dLow));
                ltl = L'*L;
                S_ki = eye(dLow) - ltl / (eye(dLow) + ltl);
                X_ki = S_ki/v * (wobs' * (ySubjObs - muSubj(obsIdx)')');
                
                % reconstruct 
                wwt = w' * w;
                if rcond(wwt) > 1e-4
                    recon = muSubj' + (w / (wwt) * (wwt + v * eye(dLow)) * X_ki)';
                else
                    rcondwarnings = rcondwarnings + 1;
                    recon = muSubj' + X_ki' * w';
                end
                    
                yRecon{i} = recon;
                
%                 yRecon{i}(obsIdx) = ySubjObs;
                yReconChk{i} = muSubj' + X_ki' * w';
                

            end
            if rcondwarnings > 0
                warning('WW'' is badly conditioned for %d patches. Forcing y = Wx + mu', rcondwarnings);
            end
            assert(all(cellfun(@isclean, yRecon)));
            assert(all(cellfun(@isclean, yReconChk)));
            
            for i = 1:N
                % put back nans where we're not sure.
                nanIdx = ~ydsmasksFullVoxels{i} & ~(data.rWeight{i} > misThr);
                yRecon{i}(nanIdx) = nan;
                yReconChk{i}(nanIdx) = nan;
            end
            
            varargout{1} = yRecon;
            varargout{2} = yReconChk;
            
        case 'latentSubspace'
            %error('Need to re-look over.');
            K = data.K;
            Y = data.Y;
            
            [dHigh, dLow, ~] = size(wg.params.W);
            
            [lll, expect] = wg.estep(data);
            % maxk = argmax(expect.gammank, [], 2);
            mi = argmax(expect.gammank, [], 2);
            
            % compute Xhat
            Xhat = zeros(size(Y,1), dLow, K);
            for i = 1:size(Y, 1)
                obsIdx = data.W(i, :) == 1;
                yobs = Y(i, obsIdx);

                for k = 1:K
                    v = wg.params.sigmasq(k);
                    w = wg.params.W(obsIdx, :, k);
                    L = w ./ sqrt(v);
                    
                    % ppca() uses Sherman-Morrison, probably because it's safer?
                    % S_ki = inv(L * L' + eye(dLow));
                    ltl = L'*L;
                    S_ki = eye(dLow) - ltl / (eye(dLow) + ltl);
                    X_ki = S_ki/v * (w' * (yobs - wg.params.mu(k, obsIdx))');
                    Xhat(i,:,k) = X_ki';
                    
                    1;
                end
            end
            assert(isclean(Xhat) && ~all(Xhat(:)==0), 'Xhat is not clean :(');
            
            yRecon = Y*0;
            xReconChk = Y*0;

            
            for k = 1:K
                w = wg.params.W(:,:,k);
                
                wwt = w'*w;
                cRecon = w / (wwt) * (wwt + wg.params.sigmasq(k) * eye(dLow)) * Xhat(mi==k, :, k)';
                yRecon(mi==k, :) = bsxfun(@plus, cRecon', wg.params.mu(k, :));
                
                xReconChk(mi==k, :) = bsxfun(@plus,  Xhat(mi==k, :, k) * w', wg.params.mu(k, :));
            end
            varargout{1} = yRecon;
            varargout{2} = xReconChk;
            
        case 'wLatentSubspace'
            %error('Need to re-look over.');
            K = data.K;
            Y = data.Y;
            
            [dHigh, dLow, ~] = size(wg.params.W);
            
            [lll, expect] = wg.estep(data);
            % maxk = argmax(expect.gammank, [], 2);
            mi = argmax(expect.gammank, [], 2);
            
            % compute Xhat
            Xhat = zeros(size(Y,1), dLow, K);
            for i = 1:size(Y, 1)
                obsIdx = data.W(i, :) == 1;
                yobs = Y(i, obsIdx);
                wtobs = data.wts(i, obsIdx);
                if sum(wtobs) == 0, 
                    continue
                end

                for k = 1:K
                    v = wg.params.sigmasq(k);
                    w = wg.params.W(obsIdx, :, k);
                    L = w ./ sqrt(v);
                    
                    % ppca() uses Sherman-Morrison, probably because it's safer?
                    % S_ki = inv(L * L' + eye(dLow));
                    ltl = L'*L;
                    S_ki = eye(dLow) - ltl / (eye(dLow) + ltl);
                    X_ki = S_ki/v * (w' * (wtobs .* (yobs - wg.params.mu(k, obsIdx)))') ./ mean(wtobs);

                    % update large matrices
                    Xhat(i,:,k) = X_ki';
                end
            end
            
            yRecon = Y*0;
            xReconChk = Y*0;

            % recon y
            for k = 1:K
                w = wg.params.W(:,:,k);
                
                wwt = w'*w;
                cRecon = w / (wwt) * (wwt + wg.params.sigmasq(k) * eye(dLow)) * Xhat(mi==k, :, k)';
                yRecon(mi==k, :) = bsxfun(@plus, cRecon', wg.params.mu(k, :));
                
                xReconChk(mi==k, :) = bsxfun(@plus,  Xhat(mi==k, :, k) * w', wg.params.mu(k, :));
            end
            varargout{1} = yRecon;
            varargout{2} = xReconChk;
            
            
        case 'latentMissingR'
            
            % data
            R = data.R;                     % R rotations from atlas to subject space.
            Y = data.Y;                     % data in subject space
            ydsmasks = data.ydsmasks;       % down-sample mask cell -- (should be "planes" in subject space)
            ydsmasksFullVoxels = data.ydsmasksFullVoxels;
            N = numel(Y);
            assert(data.K == size(wg.params.mu,1));
            
            % get maximum cluster assignment via e-step. Make sure the Estep is done in LMR
            name = wg.opts.model.name;
            wg.opts.model.name = 'latentMissingR';
            [lll, expect] = wg.estep(data);
            maxk = argmax(expect.gammank, [], 2);
            wg.opts.model.name = name;
            
            % indexing in columns is *much* faster (esp. for sparse matrices) than indexing in rows -- so we
            % transpose the R data matrix and index in columns instead of rows.
            RdataTrans = R.data';
            
            % Prepare existing sigma as a K-by-1 cell. Indexing into a small cell is faster than
            % indexing into a 3-D matrix [nData-by-nData-by-K], and since we access it this K times
            % per iteration, it ends up mattering for our iterations.
            sigmaPrevCell = dimsplit(3, wg.params.sigma);
            
            % throw a warning reminding us of a bunch of choices we're making about how to deal with data on
            % the edges. First, we only re-estimate data in the subject space that we can contribute to with
            % more than misThr weight of the atlas space voxels. Then we only use those voxels (and the
            % observed voxels) in computing Yhat (the re-estimated Y in atlas space).
            warning('Choice of which voxels to include in rotation are uncertain.');
            misThr = 0.97;
            % obsThr = 0.97; % should use for obsIdx?    
            
            rrerr = [];
            yRecon = cell(1, N);
            for i = 1:N
                k = maxk(i);
                
                % extract the observed y values in original space.
                obsIdx = ydsmasksFullVoxels{i};
%                 obsIdx = ydsmasks{i};
                misIdx = data.rWeight{i} > misThr & ~ydsmasks{i}; % used to be just ~obsIdx
%                 misIdx = ~obsIdx;
                ySubj = Y{i};
                ySubjObs = ySubj(obsIdx);
                
                % prepare the rotation matrices
                r = RdataTrans(:, R.idx{i})';
                robs = r(obsIdx, :);
                rmis = r(misIdx, :);
        
                % rotate statistics to subject space. see paffine.atl2SubjGauss().
                sigmaPrev = sigmaPrevCell{k};
                srobs = (sigmaPrev * robs');
                oosigmak = robs * srobs;
                mosigmak = rmis * srobs;

                mu = wg.params.mu(k, :);
                muSubj = mu * r';
                
                % update Y
                ySubjk = ySubj(:)';
%                 ySubjk = nan(size(ySubj(:)'));
                ySubjk(misIdx) = muSubj(misIdx) + (mosigmak * (oosigmak \ (ySubjObs - muSubj(obsIdx))'))';
                ySubjk(ydsmasks{i}) = ySubj(ydsmasks{i});
                
                % put back into output
                yRecon{i} = ySubjk;
                yRecon{i}(~misIdx & ~ydsmasks{i}) = nan;
                
%                 if rand < 0.05
%                     disp('hi');
%                 end
                
                
                
                % test recon from recon.
%                 obsIdx = ydsmasksFullVoxels{i};
%                 
%                 q = data.yrotmasks{i};
%                 q(1:round(size(q, 1)/2), :, :) = false;
%                 q = q(data.yrotmasks{i})';
%                 
% %                 f = find(obsIdx);
% %                 f = f(1:round(numel(f)/4));
% %                 q = false(size(obsIdx));
% %                 q(f) = true;
%                 % q = isodd(1:numel(obsIdx));
%                 misIdx = obsIdx & q;
%                 obsIdx = obsIdx & ~q;
%                 
%                 ySubjObs = yRecon{i}(obsIdx);
%                 
%                 % prepare the rotation matrices
%                 r = RdataTrans(:, R.idx{i})';
%                 robs = r(obsIdx, :);
%                 rmis = r(misIdx, :);
%                 sigmaPrev = sigmaPrevCell{k};
%                 srobs = (sigmaPrev * robs');
%                 oosigmak = robs * srobs;
%                 mosigmak = rmis * srobs;
%                 mmsigmak = rmis * (sigmaPrev * rmis');
%                 muSubj = mu * r';
%                 
%                 % reconstruct "original" data
%                 q = muSubj(misIdx) + (mosigmak * (oosigmak \ (ySubjObs - muSubj(obsIdx))'))';
% %                 yRecon{i} = ySubj * nan;
% %                 yRecon{i}(misIdx) = q;
%                 % compute error compared to "original" data
%                 rrerr = [rrerr; msd(q(:), ySubj(misIdx(:))')]; % reconstructed values
% %                 z = q - muSubj(misIdx);
% %                 rrerr = [rrerr; z(:)' * (mmsigmak \ z(:))];
                
            end
            varargout{1} = yRecon;
            varargout{2} = expect.gammank;
            fprintf('LogLike: %3.2f, rrerr: %3.2f, medrerr: %3.4f\n', lll, sum(rrerr), median(rrerr));
            
%             z = double(misIdx);
%             z(misIdx) = q;
%             view3Dopt({ySubj, z}, 'voxMask', data.yrotmasks{i});
%             figure(); hist(rrerr, 20); %, linspace(0, 0.3, 20));
            
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

function rpatches = reconpca(wgmm, data, prctile)
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
