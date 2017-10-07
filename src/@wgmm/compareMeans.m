function compareMeans(gmms, patchSize, titles)
    narginchk(2, 3);
    gmms = ifelse(iscell(gmms), gmms, {gmms});
    mus = cellfunc(@(x) x.params.mu, gmms);

    T = numel(gmms);
    if nargin == 1
        titles = repmat({''}, [1, numel(gmms)]);
    end
    
    % mus min and max
    mumin = min(cellfun(@(x) min(x(:)), mus));
    mumax = max(cellfun(@(x) max(x(:)), mus));
    murange = [mumin, mumax];

    % to have some sort of order, we sort each mu list with respect to the first mu in the first gmm
    mainmu = mus{1}(1, :);
    mudistord = cell(1, T);
    for i = 1:T
        mupdist = pdist2(mainmu, mus{i});
        [~, mudistord{i}] = sort(mupdist, 'ascend');
    end
    mudistord{:}
    
    % visualize mus via plots
    % [nRows, nCols] = subgrid(T);
    nRows = 1;
    nCols = T;
    figuresc();
    for i = 1:T
        subplot(nRows, nCols, i); 
        plot(mus{i}'); 
        ylim(murange); 
        axislabels('voxels', 'intensities', titles{i});
    end
    
    % visualize mus via flat images
    figuresc();
    for i = 1:T
        subplot(nRows, nCols, i); 
        imagesc(mus{i}(mudistord{i}, :), murange); 
        colormap gray; set(gca,'xticklabel',[]);
        axislabels('voxels', 'mus', titles{i});
    end

    % visualize mus via 2d slices images
    r32 = @(x) patchview.reshapeto2D(x, patchSize);
    figuresc();
    for i = 1:T
        subplot(nRows, nCols, i); 
        imagesc(r32(mus{i}(mudistord{i}, :)), murange); 
        colormap gray; axis off
        axislabels('voxels', 'mus', titles{i});
    end