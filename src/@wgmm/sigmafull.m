function [sigma, sigmainv, sigmacore, sigmarecon, sigmamerge] = ...
    sigmafull(mu, X, W, K, gammank, methods, opts, wg)
% compute the full sigma, including the sigma core (sigmac), reconstruction (sigmar), ...
% regularization and PD-fix

    % compute the core
    sigmacore = wgmm.sigmacore(mu, X, W, K, gammank, methods.core, opts, wg);
    sumgammank = sum(gammank);
    
    wtw = W' * W;
    
    sigma = zeros(size(sigmacore));
    sigmarecon = zeros(size(sigmacore));
    sigmamerge = zeros(size(sigmacore));
    sigmainv = zeros(size(sigmacore));
    for k = 1:K
        
        % compute reconstruction and merge if required (model3)
        if strcmp(methods.core, 'model3')
            % reconstruction
            sigmarecon(:,:,k) = wgmm.sigmarecon(sigmacore(:,:,k), wtw, methods.recon);
            
            % merge core with reconstruction
            switch methods.merge
                
                case 'wfact'
                    margs = {opts.mergeargs};
                case 'wfact-mult'
                    margs = {opts.mergeargs};
                case 'wfact-mult-adapt'
                    margs = {opts.mergeargs{:}, size(X, 1), sum(gammank(:, k))};
                case 'freq-prior'                    
                    margs = {opts.mergeargs{:}};
                case 'none'
                    margs = {};
                otherwise
                    error('wgmm.sigmafull: Unknown combo method');
            end
            
%             
%             if strcmp(methods.merge, 'wfact-mult-adapt'), 
%                 margs = {opts.mergeargs{:}, size(X, 1), sum(gammank(:, k))};
%             else
%                 margs = {opts.mergeargs};
%             end
            sigmamerge(:,:,k) = wgmm.sigmamerge(sigmacore(:,:,k), sigmarecon(:,:,k), ...
                wtw, methods.merge, margs{:});
            sigma(:,:,k) = sigmamerge(:,:,k);
        else
            sigma(:,:,k) = sigmacore(:,:,k);
        end

        % check PD-ness. 
        % Note: forcing pd-ness with this method seems to screw up inverses in some cases when they
        % weren't problematic originally. :(
        [~, err] = cholcov(sigma(:,:,k));
        if err > 0 || isnan(err)
            sys.warn('sigma not PD. Applying nearestSPD', 'singleWarn', true);
            sigma(:, :, k) = nearestSPD(sigma(:, :, k));
        end  

        % add diagonal regularization
        sigma(:,:,k) = sigma(:,:,k) + eye(size(sigma, 1)) * opts.sigmareg; 
        [~, err] = cholcov(sigma(:,:,k));
        assert(err == 0, sprintf('found %d/%d negative eigenvalues', err, size(sigma, 1)));

        % inverse
        [sigmainv(:,:,k), maxd] = invertsigma(sigma(:,:,k));
        if maxd > 1e-3
            warning('bad %d sigma: max(|S*invS - I|) == %3.5f, sumgammank(k)=%3.2f', ...
                k, maxd, sumgammank(k));
        end
    end
    
    % some more final checks
    assert(isclean(sigma));
    assert(isclean(sigmainv));
end

    
function [sigmainv, maxd] = invertsigma(sigma)
    sigmainv = inv(sigma);
    diff = abs(sigma * sigmainv - eye(size(sigma, 1)));
    maxd = max(diff(:)); 
end
