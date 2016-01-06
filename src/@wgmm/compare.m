function compare(gmms, patchSize, gmmMethods, varargin)
% compare the resulting gmms (and optionally reconstructions
% (gmms, patchSize, gmmMethods, varargin)
    
    inputs = parseInputs(gmms, patchSize, gmmMethods, varargin{:});
    
    % compare mus
    if ismember('mu', inputs.comparisons)
        mutitles = cellfunc(@(x) [x, '-mu'], gmmMethods);
        wgmm.compareMeans(gmms, patchSize, mutitles);
    end
    
    % compare sigmas
    if ismember('sigmas', inputs.comparisons)
        sigmatitles = cellfunc(@(x) [x, '-sigmas'], gmmMethods);
        wgmm.compareSigmas(gmms, sigmatitles);
    end
    
    % compare a single cluster
    % assumes first gmm is iso-gmm and second is iso-wgmm
    if ismember('gmm_vs_wgmm1', inputs.comparisons)
        assert(all(gmms{1}.X(:) == gmms{2}.X(:)));
        assert(~isempty(inputs.locpad));
        
        % start: get assignment from gmm, then compare sigmas.
        [~, qi] = max(gmms{1}.logpost, [], 2);
        h = hist(qi, 1:gmms{1}.K);
        [~, hi] = sort(h, 'ascend');
        
        allsigmas = {};
        alltitles = {};
        for xi = [hi(1), hi(floor(gmms{1}.K/2)), hi(gmms{1}.K)]
            xsel = gmms{1}.X(qi == xi, :);
            X = bsxfun(@minus, xsel, mean(xsel));
            Ws = gmms{2}.W(qi == xi, :);

            % compute sigmas
            fact = prod(inputs.locpad*2+1) * 15 * sum(qi == xi) / 22572;
            s1 = X' * X ./ size(X, 1);
            [s, sc, sr, sm] = model3sigma(X, Ws, 'greedy1', 'wfact-mult', fact);
            sigmas = {s1, s, sc, sr, sm};
            sigmasmall = min(cellfun(@(x) prctile(x(:), 3), sigmas));
            sigmahigh = max(cellfun(@(x) prctile(x(:), 97), sigmas));
            sigmarange = [sigmasmall, sigmahigh];        

            % compute wtw
            wtw = Ws' * Ws;
            wtw = wtw ./ max(wtw(:)) * (sigmarange(2) - sigmarange(1));
            wtw = wtw + sigmarange(1);

            % titles
            titles = {'direct sigma', 'W'' * W rescaled', 'm3sigma', 'm3s-core', 'm3s-recon', 'm3s-merge'};
            titles = cellfunc(@(x) sprintf('%s patchcount:%d', x, sum(qi == xi)), titles);
            
            % gather
            alltitles = [alltitles, titles{:}];
            allsigmas = [allsigmas, {s1, wtw, s, sc, sr, sm}];
        end
        
        % visualize
        view2D(allsigmas, 'caxis', sigmarange, ...
            'titles', alltitles, 'subgrid', [3, 6]); colormap gray;
        figuresc(); bar(1:gmms{1}.K, h);
        
        % gmm results but with W. Let's test how logpost performs.
        % to really compare this, though,w e would need a ds-wgmm.
%         wg = wgmm(gmms{1}.X, gmms{2}.W, gmms{1}.K, gmms{1}.mu, gmms{1}.sigma, gmms{1}.pi);
%         [~, qi2] = max(wg.logpost, [], 2); % assuming wgmm
%         figuresc(); scatter(qi, qi2);
    end
end



function inputs = parseInputs(gmms, patchSize, gmmMethods, varargin)
    gmms = ifelse(iscell(gmms), gmms, {gmms});
    gmmMethods = ifelse(iscell(gmmMethods), gmmMethods, {gmmMethods});

    T = numel(gmms);

    p = inputParser();
    p.addRequired('gmms', @iscell);
    p.addRequired('patchSize', @isnumeric);
    p.addRequired('gmmMethods', @iscellstr);
    p.addParameter('recons', cell(T, 1), @iscell);
    p.addParameter('reconTitles', {}, @iscell);
    p.addParameter('comparisons', {'mu', 'sigmas', 'gmm_vs_wgmm1'});
    p.addParameter('locpad', [], @isnumeric);
    p.parse(gmms, patchSize, gmmMethods, varargin{:});
    inputs = p.Results;

    % check recons and reconTitles sizes
    % recons should be T x nReconTypes
    % reconTitles should be 1 x nReconTypes
    sz = size(inputs.recons);
    assert(sz(1) == T, 'recons must be nGMMs x nReconTypes');
    if isempty(inputs.reconTitles)
        inputs.reconTitle = repmat({''}, [1, sz(2)]);
    else
        assert(sz(2) == size(inputs.reconTitles, 2), 'reconTitles must be 1 x nReconTypes');
    end
end
    