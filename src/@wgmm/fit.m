function fwgmm = fit(X, W, K, varargin)
% meant to be called from wgmmfit.

    % parse inputs
    [X, W, K, opt] = parseInputs(X, W, K, varargin{:});
    
    % go through different replicates. at each replicate, create a new wgmm object. Return the one
    % with the best log likelihood.
    
    wgmms = cell(1, opt.replicates);
    logliks = zeros(1, opt.replicates);
    for r = 1:opt.replicates

        % initialize wgmm settings
        wg = wgmm();
        if ~isempty(opt.updateMethod)
            wg.covarUpdateMethod = opt.updateMethod;
            wg.muUpdateMethod = opt.updateMethod;
            wg.logpUpdateMethod = opt.updateMethod;
        end
        wg.sigmareg = opt.regularizationValue;
        wg.sigmaopt = opt.regularizationWeight;

        % initialize fit
        wg.init(X, W, K, opt.initmethodArgs{:});
        
        % initialize counters
        llpchange = 1;
        ll = [];
        
        % First E step and ll
        [gammank, ll(1)] = wg.estep(X, W, K);
        
        % print inital ll
        printiter(wg, opt, ll, r);
        
        % iterate E-M
        ct = 1;
        while (llpchange > opt.TolFun) && ct <= opt.maxIter

            % M step
            wg.mstep(X, W, K, gammank);
            
            % E step
            [gammank, ll(ct+1)] = wg.estep(X, W, K);

            % check log likelihood
            sys.warnif(~(ll(ct+1) >= ll(ct)), sprintf('log lik went down in repl:%d iter:%d', r, ct));
            llpchange = (ll(ct+1) - ll(ct)) ./ abs(ll(ct));

            % update
            printiter(wg, opt, ll, r);
            ct = ct + 1;
        end
            
        % save update
        wgmms{r} = wg;
        logliks(r) = ll(end);
    end
    
    if opt.verbose > 1, 
        leg = arrayfunc(@(x) sprintf('replicate %d', x), 1:r);
        legend(leg);
    end
    
    % get the best replicate
    [~, mi] = max(logliks);
    fwgmm = wgmms{mi};
end

function printiter(wg, opt, ll, r)
% print information at an iteration

    ct = numel(ll) - 1;

    if ct == 0
         if opt.verbose > 0
             fprintf('\nreplicate %d\n', r);
             fprintf('%10s\t%10s\t%10s\n', 'iter', 'logp', 'percent change');
             fprintf('%10d\t%10f\t%10f\n', 0, ll(1), 1);
         end        
    else
        llpchange = (ll(ct+1) - ll(ct)) ./ abs(ll(ct));


        if opt.verbose > 0
            fprintf('%10d\t%10f\t%10f\n', ct, ll(ct+1), llpchange);

            if opt.verbose > 1
                subplot(121); cla;
                plot(ll); hold on;
                title('loglik');
                
                subplot(122); cla; 
                plot(X(:, 1), X(:, 2), '.'); 
                hold on; 
                scatter(wg.mu(:, 1), wg.mu(:, 2), '*');
                title('clustering animation of first 2 dimensions')
                drawnow;
            end
        end
    end
end

function [X, W, k, opt] = parseInputs(X, W, k, varargin)
    % parse inputs
    p = inputParser();
    p.StructExpand = true;
    
    p.addRequired('X', @ismatrix);
    p.addRequired('W', @(w) ismatrix(w) && all(size(X) == size(w)));
    p.addRequired('K', @isscalar);    
    
    p.addParameter('maxIter', 10, @isscalar);
    p.addParameter('replicates', 10, @isscalar);
    p.addParameter('initmethod', 'exemplar', @ischar); % see wgmminit(...) for options
    p.addParameter('initmethodArgs', {}, @iscell);
    p.addParameter('updateMethod', '', @ischar);
    
    p.addParameter('regularizationValue', 1e-7, @isscalar)
    p.addParameter('regularizationWeight', nan, @isscalar)
    p.addParameter('TolFun', 0.01, @(x) isscalar(x) && x <= 1 && x >= 0);
    
    p.addParameter('verbose', 2, @islogical);
    p.addParameter('debug', false, @islogical);
    
    p.parse(X, W, k, varargin{:});
    opt = p.Results;
end
