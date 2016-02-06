function [sampleX, ks] = sample(gmm, N, X, cidx)
% sample from the wgmm
%   gmm.sample(N) % N is scalar 
%   gmm.sample(N, X, cidx)
%   gmm.sample(N, X, cidx, patchSize); % TODO?
%
% optional arguments: X to compare to.


    % prep N appropriately
    narginchk(2, 5);
    K = size(gmm.mu, 1);

    % sample
    ks = sort(discretesample(gmm.pi, N), 'ascend');
    sampleXc = cell(K, 1);
    for k = 1:K
    %     w = W(i, :);
    %     sg = (1./w' * (1./w)) .* fwgmm.sigma(:, :, k);
        sampleXc{k} = mvnrnd(gmm.mu(k, :)', gmm.sigma(:, :, k), sum(ks == k));
    end
    sampleX = cat(1, sampleXc{:}); % note: only ok if ks is sorted
    
    % visualize if passed in X.
    if nargin > 2
        % show all samples as an images
        f1 = figure(); 
        for k = 1:K
            clstr = sprintf('mixture %d', k);
            figure(f1);
            subplot(K, 2, (k-1)*2 + 1); 
            imagesc(X(cidx == k, :)); colormap gray;
            title(['original data ', clstr]);
            
            subplot(K, 2, (k-1)*2 + 2); 
            imagesc(sampleXc{k}); colormap gray;
            title(['new samples ', clstr]);
        end
            
        % show first two dimensions on a scatter plot
        figure(); 
        for k = 1:K
            plot(sampleXc{k}(:, 1), sampleXc{k}(:, 2), '.'); 
            hold on; 
            plot(X(cidx == k, 1), X(cidx == k, 2), '+');
            title('first two dimensions');
        end
        legend({'samples', 'real data'});

        % show patches if given patchSize. TODO: maybe move this to subspacetools?
        if nargin > 4
           
        end 
    end
    