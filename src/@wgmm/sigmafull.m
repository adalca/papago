function [sigmaout, sigmainv] = sigmafull(wg, sigma)
% compute the full sigma, including the sigma core (sigmac), reconstruction (sigmar), ...
% regularization and PD-fix

    % prep some useful variables
    K = size(wg.expect.gammank, 2);
    
    % initialize various sigmas.
    sigmaout = sigma;
    sigmainv = zeros(size(sigma));
    for k = 1:K

        % check PD-ness. 
        % Note: forcing pd-ness with this method seems to screw up inverses in some cases when they
        % weren't problematic originally. :(
        [~, err] = cholcov(sigmaout(:,:,k));
        if err > 0 || isnan(err)
            sys.warn('sigma not PD. Applying nearestSPD', 'singleWarn', true);
            sigmaout(:, :, k) = nearestSPD(sigmaout(:, :, k));
        end  

        % add diagonal regularization
        sigmaout(:,:,k) = sigmaout(:,:,k) + eye(size(sigmaout, 1)) * wg.opts.regularizationValue; 
        [~, err] = cholcov(sigmaout(:,:,k));
        assert(err == 0, sprintf('found %d/%d negative eigenvalues', err, size(sigmaout, 1)));

        % inverse
        [sigmainv(:,:,k), maxd] = invertsigma(sigmaout(:,:,k));
        if maxd > 1e-3
            warning('bad %d sigma: max(|S*invS - I|) == %3.5f, sumgammank(k)=%3.2f', ...
                k, maxd, wg.expect.gammank(:, k));
        end
    end
    
    % some more final checks
    assert(isclean(sigmaout));
    assert(isclean(sigmainv));
end
    
function [sigmainv, maxd] = invertsigma(sigma)
    sigmainv = inv(sigma);
    diff = abs(sigma * sigmainv - eye(size(sigma, 1)));
    maxd = max(diff(:)); 
end
