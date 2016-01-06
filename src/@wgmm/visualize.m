function visualize(wgmm, varargin)
% with argument: compare. without: just visualize

    % prepare inputs
    narginchk(1, 2);
    [N, D] = size(wgmm.X);
    K = wgmm.K;

    % visualize first two dimensions
    hdim = figure(); 
    scatter(wgmm.X(:, 1), wgmm.X(:, 2)); 
    hold on; 
    scatter(wgmm.mu(:, 1), wgmm.mu(:, 2), '*'); 
    title('first two dimentions of data');

    % todo: visualize the most 
    hpcadim = figure(); 
    [coeff, score, ~, ~, explained, mu] = pca(wgmm.X);
    scatter(score(:, 1), score(:, 2)); 
    hold on;
    pmu = pcax.project(bsxfun(@minus, wgmm.mu', mu'), coeff)'; % project mus.
    scatter(pmu(:, 1), pmu(:, 2), '*'); 
    title(sprintf('first two (overall) pca dims (%3.2f%% varexplained)', sum(explained(1:2))));
    
    % centers
    hcen = figure(); 
    plot(wgmm.mu', 'o-'); 
    ylabel('feature / dimension');
    title('mus');
    
    % covariances
    hcov = figure(); 
    imagesc(reshape(wgmm.sigma, [D*K, D])); 
    colormap gray; axis equal off;
    title('stacked sigmas');
    
    % plot original
    if nargin == 2
        ogmm = varargin{1};
        
        % first dimensions
        figure(hdim); hold on;
        scatter(ogmm.mu(:, 1), ogmm.mu(:, 2), '+'); 

        % first pca dimensions
        figure(hpcadim); hold on;
        pmu = pcax.project(bsxfun(@minus, ogmm.mu', mu'), coeff)'; % project mus.
        scatter(pmu(:, 1), pmu(:, 2), '+'); 
        
        % centers
        figure(hcen); hold on;
        plot(ogmm.mu', '+-'), hold on;
        legend({'discovered mu', 'given mu'});
        
        % covariances. 
        [~, cequiv] = min(pdist2(wgmm.mu, ogmm.mu));
        figure(hcov); clf; 
        subplot(1, 2, 1); 
        imagesc(reshape(wgmm.sigma(:,:,cequiv), [D*K, D])); colorbar;
        colormap gray; axis equal off;
        title('stacked discovered sigmas');
        subplot(1, 2, 2);
        imagesc(reshape(ogmm.sigma, [D*K, D])); colorbar;
        colormap gray; axis equal off;    
        title('stacked given sigmas');

    end
    