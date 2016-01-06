function sigma = sigmacore(mu, X, W, K, gammank, coremethod, coreargs)

    Nk = sum(gammank, 1);
    nDims = size(X, 2);

    % sigma update
    sigma = zeros([nDims, nDims]);
    
    switch coremethod
        
        case 'model0' % no weights, standard sigma
            % compute sigma without weights
            sigma = zeros(nDims, nDims, K);
            
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :));
                
                % Variant 1 -- this is slower
                % q = 0;
                % for i = 1:size(gammank, 1)
                %     q = q + gammank(i, k)' * (xc(:, i)' * xc(:, i));
                % end
                % sigma(:,:,k) = 1 ./ sum(gammank(:, k)) * q;
                
                % Variant 2 -- this is significantly faster!
                gw = sqrt(gammank(:, k));
                wx = bsxfun(@times, gw, xc);
                aa = wx' * wx;
                sigma(:, :, k) = aa ./ sum(gammank(:, k));
            end
        
        case 'model1' 
            % model 1 in terms of derivations during the wgmm project. This implementation is using
            % matrix algebra, but running into memory issues, and not that much faster. see
            % memsafe-model1 method. This version uses Paolo de Leva's multiprod() function.
            iwAiw = coreargs.iwAiw;
            
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :))';

                 % method 1. 
                diff = permute(xc, [1, 3, 2]);
                xxt = multiprod(diff, permute(diff, [2, 1, 3])); % maybe just loop?
                ixxti = iwAiw(W, xxt); 
                s1 = sum(bsxfun(@times, ixxti, permute(gammank(:, k), [2, 3, 1])), 3);
                sigma = 1 ./ Nk(k) .* s1;
            end
                
        case 'memsafe-model1' 
            % result should be the same as 'model1'
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :))';

                % method 2. loop.
                numer = 0; 
                for i = 1:size(X, 1)
                    w = W(i, :)';
                    wx = w .* xc(:, i); 
                    q = wx * wx'; 
                    numer = numer + gammank(i, k) * q;
                end
                sigma(:,:,k) = numer ./ Nk(k);
            end
            
        case 'model3'
            % model 3, as computed in wgmm timeline. using matrix algebra, which seems significant
            % faster than the memsage method. It doesn't use too much memory, since we don't have to
            % compute as large matrices.
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :));
                
                % Variant 1 -- this is slower
                % wx = W .* xc;
                % aa = bsxfun(@times, gammank(:, k), wx)' * wx;
                % wtsum = bsxfun(@times, gammank(:, k), W)' * W;
                % sigma(:, :, k) = aa ./ wtsum;
                
                % Variant 2 -- this is significantly faster!
                gw = bsxfun(@times, sqrt(gammank(:, k)), W);
                wx = gw .* xc;
                aa = wx' * wx;
                wtsum = gw' * gw;
                sigma(:, :, k) = aa ./ wtsum;
            end
            
        case 'memsafe-model3'
            % memory safe (looping) model 3
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :))';
                numer = 0;
                denom = 0;
                for i = 1:size(X, 1)
                    w = sqrt(W(i, :)');
                    z = w .* xc(:, i);
                    q = z * z';
                    numer = numer + gammank(i, k) * q;
                    
                    denom = denom + gammank(i, k) .* (w * w');
                end
                sigma(:, :, k) = numer ./ denom;
            end
            
        otherwise
            error('unknown covarianceupdate');
    end
    
    assert(isclean(sigma));
end
