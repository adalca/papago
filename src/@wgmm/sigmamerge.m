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
            
        case 'none'
            sigma = sigmac;
            
        otherwise
            error('wgmm.sigmafull: Unknown combo method');
    end
    
    assert(isclean(sigma));