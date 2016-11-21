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
        %if ~isfield(wg.expect, 'gammank')
        if opts.maxIter > 0
            [ll(1), wg.expect] = wg.estep(data);
        else 
            ll = 0;
           %ll(1) = wg.logp(data);
        end
        wg.stats(1).expect = wg.expect;
        
        % print inital ll
        printiter(wg, opts, ll, r);
        
        % Iterative E.M. updates
        ct = 1;
        while (llpchange > opts.TolFun) && ct <= opts.maxIter
            tic;

            % M step
            wg.params = wg.mstep(data);
            wg.stats(ct+1).params = wg.params;

            % E step
            [ll(ct+1), expect, wg] = wg.estep(data);
            wg.expect = expect;
            wg.stats(ct+1).expect = wg.expect;
            
            % add time stats
            wg.stats(ct+1).toc = toc;
            wg.stats(ct+1).ll = ll(ct+1);

            % check log likelihood
            sys.warnif(~(ll(ct+1) >= ll(ct)), sprintf('log lik went down in repl:%d iter:%d', r, ct));
            llpchange = abs(ll(ct+1) - ll(ct)) ./ abs(ll(ct));

            % update
            printiter(wg, opts, ll, r);
            ct = ct + 1;
            
            % temp save.
            % save(tempname, 'wg', '-v7.3');
        end
            
        % save update
        wgmms{r} = wg;
        logliks(r) = ll(end);
    end
    
    if opts.verbose > 1
        leg = arrayfunc(@(x) sprintf('replicate %d', x), 1:r);
        legend(leg);
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
    p.addParameter('regularizationValue', dopts.regularizationValue, @isscalar)    
    p.addParameter('TolFun', dopts.TolFun, @(x) isscalar(x) && x <= inf && x >= -inf);
    p.addParameter('verbose', dopts.verbose, @isIntegerValue);
    
    % parse
    p.parse(data, varargin{:});
    
    % copy
    opts = p.Results;
    opts = rmfield(opts, 'modelName');
    opts = rmfield(opts, 'modelArgs');
    opts.model = p.Results.modelArgs;
    opts.model.name = p.Results.modelName;
    opts = rmfield(opts, 'init');
    opts = rmfield(opts, 'initArgs');
    opts.init = p.Results.initArgs;
    opts.init.method = p.Results.init;
end

function printiter(wg, opts, ll, r)
% print information at an iteration

    ct = numel(ll) - 1;

    if ct == 0
         if opts.verbose > 0
             fprintf('\nreplicate %d\n', r);
             fprintf('%10s\t%10s\t%10s\n', 'iter', 'logp', 'percent change');
             fprintf('%10d\t%10f\t%10f\n', 0, ll(1), 1);
         end
         
    else
        llpchange = (ll(ct+1) - ll(ct)) ./ abs(ll(ct));


        if opts.verbose > 0
            fprintf('%10d\t%10f\t%10f\n', ct, ll(ct+1), llpchange);

%             if opt.verbose > 1
%                 subplot(121); cla;
%                 plot(ll); hold on;
%                 title('loglik');
%                 
%                 subplot(122); cla; 
%                 plot(X(:, 1), X(:, 2), '.'); 
%                 hold on; 
%                 scatter(wg.mu(:, 1), wg.mu(:, 2), '*');
%                 title('clustering animation of first 2 dimensions')
%                 drawnow;
%             end
        end
    end
end

% note: while transfering to new code, need to remember to deal with the following setting transfers
% add defaults
% initmethod = 'exemplar';
% debug = false;
% sigmareg = nan;
% sigmaopt = {0};


% covarUpdateMethod = 'model5';
% muUpdateMethod = 'model5';
% logpUpdateMethod = 'model5';
% model4fn = @(w) diag((-log(w)).^ 2);
% 
% covarReconMethod = 'greedy1';
% covarMergeMethod = 'wfact-mult-adapt'; % 'none'
% 
% if ~isempty(opt.updateMethod)
% wg.covarUpdateMethod = opt.updateMethod;
% wg.muUpdateMethod = opt.updateMethod;
% wg.logpUpdateMethod = opt.updateMethod;
% end
% 
% p.addParameter('regularizationWeight', nan, @(x) isscalar(x) | iscell(x))
% %p.addParameter('covarMergeMethod', 'wfact-mult-adapt', @ischar);
% p.addParameter('covarMergeMethod', 'none', @ischar);
% 
% p.addParameter('model4fn', @(w) diag((-log(w)).^ 2), @isfunc);
% p.addParameter('debug', false, @islogical);
