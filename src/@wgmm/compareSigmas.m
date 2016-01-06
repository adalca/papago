function compareSigmas(gmms, titles)
    narginchk(1, 2);
    gmms = ifelse(iscell(gmms), gmms, {gmms});
    sigmas = cellfunc(@(x) x.sigma, gmms);
    [S, ~, K] = size(sigmas{1});

    T = numel(gmms);
    if nargin == 1
        titles = repmat({''}, [1, numel(gmms)]);
    end
    assert(size(titles, 2) == T);
    if size(titles, 1) == 1 && K > 1
        pt = titles;
        for k = 1:K
            for t = 1:T
                titles{k, t} = sprintf('%s K:%d', pt{1, t}, k);
            end
        end
    end
    
    % sigmas min and max
    sigmasmall = min(cellfun(@(x) prctile(x(:), 3), sigmas));
    sigmahigh = max(cellfun(@(x) prctile(x(:), 97), sigmas));
    sigmarange = [sigmasmall, sigmahigh];
    
    % to have some sort of order, we sort each sigma list with respect to the first sigma in the first gmm
    mainsigma = sigmas{1}(:,:,1);
    mainsigma = mainsigma(:)';
    sigmadistord = cell(1, T);
    for c = 1:T
        r = reshape(sigmas{c}, [S*S, K])';
        sigmapdist = pdist2(mainsigma, r);
        [~, sigmadistord{c}] = sort(sigmapdist, 'ascend');
    end
    
    % visualize sigmas via flat images
    [nRl, nCols] = subgrid(K);
    nRows = nRl * T;
    
    figuresc();
    for r = 1:nRows
        for c = 1:nCols
            t = floor((r - 1) ./ nRl) + 1;
            k = mod(r - 1, nRl) * nCols + c;
            if k > K, continue; end
            
            idx = sub2ind([nCols, nRows], c, r);
            subplot(nRows, nCols, idx); 
            
            imagesc(sigmas{t}(:, :, sigmadistord{t}(k)), sigmarange); 
            colormap gray; axis off;
            title(titles{k, t});
        end
    end
    