function sigma = sigmacore(mu, X, W, K, gammank, coremethod, coreargs, wg)

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
        
        case 'model4exp'
            
            for k = 1:K
                xc = bsxfun(@minus, X, mu(k, :))';
                
                numer = 0;
                denom = 0;
                for i = 1:size(X, 1)
                    w = sqrt(W(i, :)');
                    Di = wg.model4fn(w);
                    s = wg.sigma(:,:,k);
                    sg = s + Di;
                    xtx = xc(:, i) * xc(:, i)';
                    
                    numer = numer + gammank(i, k) * (xtx + (s / sg) * Di);
                    
                    denom = denom + gammank(i, k);
                end
                sigma(:, :, k) = numer ./ denom;
                imagesc(sigma(:, :, k)); 
                drawnow;
                pause(0.01);
            end
            
            
        case 'model5'
%             for k = 1:K
%                 muk = wg.mu(k, :);
%                 sigma(:,:,k) = fminsearch(@(s) wg.model5exp(s, X, W, muk, wg.model4fn), wg.sigma(:,:,k));
%             end
            
%             for k = 1:K
%                 numer = 0;
%                 denom = 0;
%                 sigmak = wg.sigma(:,:,k);
%                 for i = 1:size(X, 1)
%                     w = W(i, :);
%                     x = X(i, :);
%                     df = (x + wg.mu(k,:));
%                     dfd = df(:) * (df(:))';
%                     
%                     % sigmas
%                     Di = wg.model4fn(w);
%                     sg = sigmak + Di;
%                     
%                     numer = numer + gammank(i, k) * dfd / (sg + Di) / Di / (inv(sigmak) + inv(Di));
%                     denom = denom + gammank(i, k) * inv(sg + Di);
%                 end
%                 sigma(:, :, k) = denom \ numer; 
%             end    
            
%             for k = 1:K
%                 numer = 0;
%                 denom = 0;
%                 sigmak = wg.sigma(:,:,k);
%                 for i = 1:size(X, 1)
%                     w = W(i, :);
%                     x = X(i, :);
%                     df = (x - wg.mu(k,:));
%                     dfd = df(:) * (df(:))';
%                     
%                     % sigmas
%                     Di = wg.model4fn(w);
%                     sg = sigmak + Di;
%                     
%                     numer = numer + gammank(i, k) * (- eye(numel(x)) + (dfd / (sg + Di))) / (sg + Di);
%                 end
%                 sigma(:, :, k) = sigmak - 0.00001 * numer ./ sum(gammank(:, k)); 
%             end    

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % MICCAI2016 Push implementation.
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fprintf(2, 'current using init with sigma of y (noisy data). Maybe use previous sigma?');
% 
%         for k = 1:K
%             muk = wg.mu(k, :);
%             xc = bsxfun(@minus, X, muk);
%             
%             % compute sigma for this cluster ffrom noisy data xc.
%             gw = sqrt(gammank(:, k));
%             wx = bsxfun(@times, gw, xc);
%             aa = wx' * wx; % gamma * (x-u) (x-u)'
%             sigmainitk = aa ./ sum(gammank(:, k));
%             
%             % hack 1
%             numer = 0; denom = 0; 
%             meanD = 0;
%             for i = 1:size(X, 1)
%                 w = W(i, :); 
%                 x = X(i, :);
%                 df = (x - muk);
%                 dfd = df(:) * (df(:))';
% 
%                 Di = wg.model4fn(w);
%                 sg = sigmainitk + Di;
% 
%                 numer = numer + (sg \ dfd);
%                 %denom = denom + inv(sg);
%                 denom = denom + inv(sg) + inv(sg) * dfd * inv(sigmainitk) * Di * inv(sigmainitk);
% 
%                 meanD = meanD + Di;
%             end
%             meanD = meanD ./ size(X, 1);
%             sigma(:,:,k) = (denom \ numer) - meanD ./ 8; % empirical ./4
%             % sigma(:,:,k) = (denom \ numer); 
%             
%         end

            for k = 1:K
                muk = wg.mu(k, :);

                % model 5
                df = bsxfun(@minus, wg.expect.Xk(:,:,k), muk);
                numer = 0; 
                for i = 1:size(X, 1)
                    % extract important terms.
                    w = W(i, :);
                    Di = wg.model4fn(w);
                    dfd = df(i, :)' * df(i, :);
                    sk = wg.sigma(:,:,k);

                    % previous sigma^{Rt}_{ik}. We can't precompute this due to the size.
                    sigmaikr = Di * ( (Di + sk) \ sk);

                    % update numerator and denominator.
                    numer = numer + gammank(i, k) .* (sigmaikr + dfd);
                end
                sigma(:,:,k) = numer ./ sum(gammank(:, k));
            end
        
        otherwise
            error('unknown covarianceupdate');
    end
    
    assert(isclean(sigma));
end
