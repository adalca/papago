function fwgmm = fit(data, varargin)
% meant to be called from wgmmfit.

    % parse inputs
    [data, opts] = parseInputs(data, varargin{:});

    mastertic = tic;
    
    % go through different replicates. at each replicate, create a new wgmm object. 
    % Return the one with the best log likelihood.
    wgmms = cell(1, opts.replicates);
    logliks = zeros(1, opts.replicates);
    for r = 1:opts.replicates

        % construct wgmm
        wg = wgmm(opts);        

        % initialize fit
        wg.init(data, opts.init);
        wg.stats(1).params = wg.params;
        wg.opts.model
        
        % initialize counters
        llpchange = 1;
        ll = [];
        
        % First E step and ll
        if opts.maxIter > 0 % only compute expectation and ll
            [ll(1), wg.expect] = wg.estep(data);
        else 
            ll = 0;
        end
        wg.stats(1).ll = ll;
        wg.stats(1).expect = wg.expect;
        wg.stats(1).params = wg.params;
        
        % print inital ll
        if wg.opts.verbose > 0, printiter(wg, ll, r); end
        
        % Iterative E.M. updates
        ct = 1;
        while (llpchange > opts.TolFun) && ct <= opts.maxIter
            itertic = tic;

            % recluster if necessary
            wg = wg.recluster();
            
            % M step
            mtic = tic;
            wg.params = wg.mstep(data);
            wg.stats(ct+1).params = wg.params;
            wg.stats(ct+1).mtoc = toc(mtic);

            % E step
            etic = tic;
            [ll(ct + 1), expect] = wg.estep(data);
            wg.expect = expect;
            wg.stats(ct+1).expect = wg.expect;
            wg.stats(ct+1).etoc = toc(etic);
            
            % add time stats
            wg.stats(ct + 1).toc = toc(itertic);
            wg.stats(ct + 1).ll = ll(ct+1);

            % check log likelihood
            if(ll(ct + 1) < ll(ct)) 
                fprintf(2, 'wgmm.fit: log likelihood went down in repl:%d iter:%d\n', r, ct);
            end
            llpchange = abs(ll(ct+1) - ll(ct)) ./ abs(ll(ct));

            % verbosity
            if wg.opts.verbose > 0, printiter(wg, ll, r); end
            
            % update
            ct = ct + 1;
            
            % temp save.
            % save(tempname, 'wg', '-v7.3');
        end
            
        % save update
        wgmms{r} = wg;
        logliks(r) = ll(end);
    end
    
    % get the best replicate
    [~, mi] = max(logliks);
    fwgmm = wgmms{mi};
    
    % add master toc
    fwgmm.stats(end).mastertoc = toc(mastertic);
end

function [data, opts] = parseInputs(data, varargin)
    % get defaults;
    dopts = wgmm.optionDefaults();

    % parse inputs
    p = inputParser();
    p.StructExpand = true;
    
    % data
%     p.addRequired('X', @ismatrix);
%     p.addRequired('W', @(w) ismatrix(w) && all(size(X) == size(w)));
%     p.addRequired('K', @isscalar);    
    p.addRequired('data', @isstruct);    
    
    % model options
    p.addParameter('modelName', dopts.model.name, @ischar);
    p.addParameter('modelArgs', dopts.model, @isstruct);
    
    % initialization
    p.addParameter('init', dopts.init.method, @ischar); % see wgmminit(...) for options
    p.addParameter('initArgs', dopts.init, @isstruct);
    
    % fitting options
    p.addParameter('maxIter', dopts.maxIter, @isscalar);
    p.addParameter('replicates', dopts.replicates, @isscalar);
    p.addParameter('reclusterThreshold', dopts.reclusterThreshold, @isscalar);
    p.addParameter('reclusterMethod', dopts.reclusterMethod, @isscalar);
    p.addParameter('regularizationValue', dopts.regularizationValue, @isscalar)    
    p.addParameter('TolFun', dopts.TolFun, @(x) isscalar(x) && x <= inf && x >= -inf);
    p.addParameter('verbose', dopts.verbose, @isIntegerValue);
    
    % parse
    p.parse(data, varargin{:});
    
    % copy
    opts = p.Results;
    opts = rmfield(opts, 'data');
    opts = rmfield(opts, 'modelName');
    opts = rmfield(opts, 'modelArgs');
    opts.model = p.Results.modelArgs;
    opts.model.name = p.Results.modelName;
    opts = rmfield(opts, 'init');
    opts = rmfield(opts, 'initArgs');
    opts.init = p.Results.initArgs;
    opts.init.method = p.Results.init;
end

function printiter(wg, ll, r)
% print information at an iteration
    
    ct = numel(ll) - 1;

    if ct == 0
         fprintf('\nreplicate %d\n', r);
         fprintf('%5s\t%15s\t%10s\t%25s\n', ...
             'iter', 'logp', 'perc. change', 'time: total | mstep | estep');
         fprintf('%5d\t%15s\t%10s\t%27s\n', 0, num2bank(ll(1)), ' ', ' ');
         
    else
        % compute change
        llpchange = (ll(ct+1) - ll(ct)) ./ abs(ll(ct));

        % compute the number of decimal digits necessary
        tolfun = wg.opts.TolFun;
        % Don't count sig digits before decimal pt
        % https://www.mathworks.com/matlabcentral/answers/
        %   /184255-how-to-calculate-digits-after-decimal-point?s_tid=answers_rc1-3_p3_Topic
        llChangeDigits = max(sigdigits(tolfun)-floor(log10(abs(tolfun)))-1,0); 
        llChangeStr = sprintf('%%9.%df%%%%', llChangeDigits+1);
        
        % print
        fmt = ['%5d\t%15s\t', llChangeStr, '\t%7.1f | %7.1f | %7.1f\n'];
        fprintf(fmt, ct, num2bank(ll(ct+1)), llpchange, ...
            wg.stats(end).toc, wg.stats(end).mtoc, wg.stats(end).etoc);
    end
end
