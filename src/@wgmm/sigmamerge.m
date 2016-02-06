function sigma = sigmamerge(sigmac, sigmar, wtw, method, varargin)
% merge sigma core with sigma reconstruction estimate via weights.

    switch method
        case 'wfact'
            fact = varargin{:}; % in experiments, used (prod(locpad * 2 + 1) * 15 ./ K)
            ww = min(wtw, fact) ./ fact;
            sigma = ww .* sigmac + (1-ww) .* sigmar;
            
        case 'wfact-mult'
            fact = varargin{:}; % in experiments, used (prod(locpad * 2 + 1) * 15 ./ K)
            ww = min(wtw, fact) ./ fact;
            mult = median(sigmac(ww == 1) ./ sigmar(ww == 1));
            assert(isclean(mult));
            sigma = ww .* sigmac + (1-ww) .* sigmar .* mult;
            
        case 'wfact-mult-adapt'
            fact = varargin{1}; % in experiments, used (prod(locpad * 2 + 1) * 15 ./ K)
            nTotalSubj = varargin{2};
            nSubj = varargin{3}; % should be nSubj/nTotalSubj
            fact = fact .* nSubj / nTotalSubj;
            assert(isclean(fact));
            ww = min(wtw, fact) ./ fact;
            mult = median(sigmac(ww == 1) ./ sigmar(ww == 1));
            if ~isclean(mult)
                warning('mult is not clean: %f', mult);
                mult = median(sigmac(ww == max(ww(:))) ./ sigmar(ww == max(ww(:))));
            end
            sigma = ww .* sigmac + (1-ww) .* sigmar .* mult;
            
        case 'freq-prior'
            fact = varargin{1};
            patchSize = varargin{2};
            maxIter = 1; % hardcoded for now. TODO maybe change?
            maxFreq = maxIter*(2*pi/ max(patchSize) );
            %maxFreq2 = maxIter*(2*pi./ patchSize );
            
            [x, y, z] = ndgrid(1:patchSize(1), 1:patchSize(2), 1:patchSize(3));
            ww = min(wtw, fact) ./ fact;
            
            sincFun = zeros(prod(patchSize), prod(patchSize));
            %sincFun2 = zeros(prod(patchSize), prod(patchSize));
            
            for k=1:prod(patchSize)
                sincFun(k,:) = sinc(maxFreq/pi * pdist2([x(k) y(k) z(k)], [x(:) y(:) z(:)]) );
                %sincFun2(k,:) = sinc(maxFreq2(1)/pi * abs(x(k)-x(:)) ) .*  ...
                %    sinc(maxFreq2(2)/pi * abs(y(k)-y(:)) ) .* ...
                %    sinc(maxFreq2(3)/pi * abs(z(k)-y(:)) ); 
            end
            sincFun = sincFun .* mean(diag(sigmac));
            sigma = ww.*sigmac + (1-ww).*sincFun;
            
        case 'none'
            sigma = sigmac;
            
        otherwise
            error('wgmm.sigmafull: Unknown combo method');
    end
    
    assert(isclean(sigma));
    