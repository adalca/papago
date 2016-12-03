function [wg, wgLS, wgGreedy, wgLSGreedy] = fitgmdist2wgmmLS(y, K, regVal, reps, gmmopt, dLow)
% fit a gmm distribution and transform it to a wgmm, and a wgmm-LS (dLow-based parameters)
%
% wg - the gmdist cast to a wgmm
% wgLS - the wgmm with W computed from wg params
% greedy - a wgmm with params computed from argmax of each cluster

    dHigh = size(y, 2);

    % fit matlab gmdist
    gmdist = fitgmdist(y, K, 'regularizationValue', regVal, 'replicates', reps, 'Options', gmmopt);
    
    % compute expectations
    dspost = gmdist.posterior(y);
    [~, dsmi] = max(dspost, [], 2);
    wg = wgmm.gmdist2wgmm(gmdist);
    wg.expect.gammank = dspost;

    % get greedy wg
    wgGreedy = wgmm(); % greedy
    for k = 1:K
        X0orig = y(dsmi == k, :);
        wgGreedy.params.mu(k, :) = mean(X0orig);
        wgGreedy.params.sigma(:,:,k) = cov(X0orig);
        wgGreedy.params.pi(k) = sum(dsmi == k) ./ numel(dsmi);
    end
    wgGreedy.params.pi = wg.params.pi;
    
    % get a wg-LS
    if exist('dLow', 'var')
        wgLS = wgmm();
        lsW = zeros(dHigh, dLow, K);
        lsv = zeros(1, K);
        dsNewSigma = zeros(dHigh, dHigh, K);
        for k = 1:K
            % compute W, v, etc
            [lsW(:,:,k), lsv(k), dsNewSigma(:,:,k)] = sigma2modelParams(wg.params.sigma(:,:,k), dLow);
            wgLS.params.sigma(:, :, k) = lsW(:,:,k) * lsW(:,:,k)' + lsv(k) * eye(dHigh);
        end
        wgLS.params.pi = wg.params.pi;
        wgLS.params.W = lsW;
        wgLS.params.sigmasq = lsv;
        wgLS.params.mu = wg.params.mu;


        % get a wg-LS
        wgLSGreedy = wgmm();
        lsW = zeros(dHigh, dLow, K);
        lsv = zeros(1, K);
        dsNewSigma = zeros(dHigh, dHigh, K);
        for k = 1:K
            % compute W, v, etc
            [lsW(:,:,k), lsv(k), dsNewSigma(:,:,k)] = sigma2modelParams(wgGreedy.params.sigma(:,:,k), dLow);
            wgLSGreedy.params.sigma(:, :, k) = lsW(:,:,k) * lsW(:,:,k)' + lsv(k) * eye(dHigh);
        end
        wgLSGreedy.params.pi = wgGreedy.params.pi;
        wgLSGreedy.params.W = lsW;
        wgLSGreedy.params.sigmasq = lsv;
        wgLSGreedy.params.mu = wgGreedy.params.mu;
    else
        wgLs = [];
        wgLsGreedy = [];
    end
