%% xplorSigmaReconstruction
% exploration of sampling of $P(x_3|x_1)$ when we don't have the covariance of
% $x_3$ and $x_1$, but have the covariance of $x_3$ -and- $x_k$, and $x_k$ -and-
% $x_1$, for some $x_k$.
%
% We will use the following terminology:
%
% * "intermediate" refers to xk, or sometimes referred to x2
%
% * sample $\sigma_{13}$ is the covariance obtained by first sampling from the process 
%   $p(x_3|x_1)=\sum_k p(x_3|x_k)p(x_k|x_1)$, and computing the covariance of the resulting
%   samples
% * data $\sigma_{13}$ is the covariance between $x_1$ and $x_3$ computed from real data 
% * formula $\sigma_{13}$ is the covariance between $x_1$ and $x_3$ computed via our formula
%   (wmean1)
% * formula v2 $\sigma_{13}$ is the covariance between $x_1$ and $x_3$ computed via our second
%   formula (greedy1) for greedy1, we basically choose a x_k that has high covariance with x_1
%   (weight appropriately, and taking into account the patch weights in real data), and using the
%   main formula using just that xk (no mean)
%
% * iso refers to data/patches from full isotropic data. 
% * biso refers to data from blurred isotropic volumes, with spatial sigma=2 (pretty blurry)
% * bisow refers to data from combined blurry + iso data. a bisow volume is w * iso + (1-w) * biso
% where w are the mask volume/patches.

%% setup

% visualization settings
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

sampleCovColorSet = {'MarkerEdgeColor', [0, 0.447, 0.741], 'MarkerEdgeColor', [0, 0.447, 0.741]};
formulaCovColorSet = {'MarkerEdgeColor', [0.85, 0.325, 0.098], 'MarkerEdgeColor', [0.85, 0.325, 0.098]};
dataCovColorSet = {'MarkerEdgeColor', [0.929, 0.694, 0.125], 'MarkerEdgeColor', [0.929, 0.694, 0.125]};
sampleCovIndivColorSet = {'MarkerEdgeColor', [0.5, 0.75, 1], 'MarkerFaceColor', [0.5, 0.75, 1]};

% realDataFile. needs to have isopatches, bisopatches, maskpatches
warning('TODO: Get new data, and diff locs. Current data is old. Also, run ds test, but just bisow');
realDataFile = [SYNTHESIS_DATA_PATH, 'adni_locpad1']; 

%% (single intermediate): changing just one covariation (element of sigma) at a time
% vary the jth element of sigma to see what happens.
%
% here, the sigma vector is 1x6 and contains 
%   $[\sigma_{kk}, \sigma_{1k}, \sigma_{11}, \sigma_{33}, \sigma_{k3}, \sigma_{kk}]$

% parameters
nSamples = 100000; 
nVarySteps = 22; % number of steps for varying sigma_j
jMax = 5;
otherSigmaMax = 10;

% prepare the figure
figuresc();
[nRows, nCols] = subgrid(6);

% go through each sigma element
for j = 1:6;
    
    % initial sigmas
    sigma = ones(1, 6) * otherSigmaMax;
    
    % span sigma_j
    if ismember(j, [1, 3, 4, 6])
        sigmajRange = linspace(0.1, jMax, nVarySteps);
    else
        sigmajRange = linspace(-jMax, jMax, nVarySteps);
    end
    
    samplecov = zeros(nVarySteps, 1);
    formulacov = zeros(nVarySteps, 1);
    for u = 1:nVarySteps
        % vary the jth entry of sigma
        sigma(j) = sigmajRange(u);

        % sample
        [~, covSampleResult] = sampleSigmaRecon(nSamples, sigma);
        samplecov(u) = covSampleResult(3);

        % guess using our formula
        formulacov(u) = sigma(2) * sigma(5) ./ sigma(6);
    end
    
    % plot
    subplot(nRows, nCols, j); grid on;
    sigmaCovars = [22, 12, 11, 33, 23, 22];
    plot(sigmajRange, samplecov, '*-', sampleCovColorSet{:}); hold on;
    plot(sigmajRange, formulacov, 'o--', formulaCovColorSet{:});
    ylim([min(samplecov)-1, max(samplecov)+1]);
    sigmastr = sprintf('$\\sigma_{%d}$', sigmaCovars(j));
    title(['Varying ', sigmastr]);
    legend({'sample $\sigma_{13}$', 'formula $\sigma_{13}$'});
    xlabel(sigmastr);
    ylabel('estimated $\sigma_{13}$');
end

%% (single intermediate): completely randomize entire sigma vector
nTimes = 100;
nSamples = 100000;
sigmaMax = 10;

samplecov = zeros(nTimes, 1);
formulacov = zeros(nTimes, 1);
for u = 1:nTimes
    sigma = rand(1, 6)*sigmaMax;
    [~, covSampleResult] = sampleSigmaRecon(nSamples, sigma);
    samplecov(u) = covSampleResult(3);

    % guess using our formula
    formulacov(u) = sigma(2) * sigma(5) ./ sigma(6);
end
[~, si] = sort(samplecov);

% visualization
figuresc(); hold on; ax = gca; ax.XMinorGrid = 'on';
plot(1:nTimes, samplecov(si), '*', sampleCovColorSet{:}); 
plot(1:nTimes, formulacov(si), 'o', formulaCovColorSet{:}); 
legend({'sample-$\sigma_{13}$', 'formula-$\sigma_{13}$'});
xlabel('runs (sorted by sample-$\sigma_{13}$)');
ylabel('estimated $\sigma_{13}$');
title('sigma estimates for different runs');

%% prepare real data

% requires isopatches
load(realDataFile);
bisowpatches = isopatches .* maskpatches + bisopatches .* (1 - maskpatches);

% compute covariances, etc Note: can't use cov() for maskSigma since cov() first substracts the mean
isopatchesCen = bsxfun(@minus, isopatches, mean(isopatches));
bisopatchesCen = bsxfun(@minus, bisopatches, mean(bisopatches));
bisowpatchesCen = bsxfun(@minus, bisowpatches, mean(bisowpatches));
dataSigma = cov(isopatchesCen);
bdataSigma = cov(bisopatchesCen);
maskSigma = maskpatches' * maskpatches;

% visualize
figuresc();
p95 = prctile([dataSigma(:); bdataSigma(:)], 95);
subplot(221); imagesc(dataSigma, [0, p95]); colormap gray; title('covariance of iso data');
subplot(222); imagesc(bdataSigma, [0, p95]); colormap gray; title('covariance of blurry data');
subplot(223); imagesc(maskSigma); colormap gray; title('mask "covariance"');
subplot(224); imagesc(cov(bisowpatchesCen), [0, p95]); colormap gray; title('covariance of b-iso-w data');

%% (single intermediate): sampling triplets real data
% plot the detected covariance of samples p(x3|x1) vs the guessed covariance using our formula

% parameters
nFeats = size(dataSigma, 1);
nTimes = 100;
nSamples = 1000;

samplecov = zeros(nTimes, 1);
formulacov = zeros(nTimes, 1);
datacov = zeros(nTimes, 1);
for i = 1:nTimes
    % these are points 1,2,3 in our notation above. 2 is the intermediate we use.
    triplet = randsample(nFeats, 3); 

    % sample
    [~, covSampleResult, ~, sigma] = sampleSigmaRecon(nSamples, dataSigma, triplet);
    
    % compute sigma estimates
    samplecov(i) = covSampleResult(3); % sampling
    formulacov(i) = sigma(2) * sigma(5) ./ sigma(6);  % formula
    datacov(i) = dataSigma(triplet(1), triplet(3)); % real data
end
[~, si] = sort(samplecov);

% visualize
figuresc(); hold on; ax = gca; ax.XMinorGrid = 'on';
plot(1:nTimes, samplecov(si), '*', sampleCovColorSet{:}); 
plot(1:nTimes, formulacov(si), 'o', formulaCovColorSet{:}); 
plot(1:nTimes, datacov(si), 'd', dataCovColorSet{:}); 
legend({'sample-$\sigma_{13}$', 'formula-$\sigma_{13}$', 'data-$\sigma_{13}$'});
xlabel('runs (sorted by sample-$\sigma_{13}$)');
ylabel('estimated $\sigma_{13}$');
title('sigma estimates for different runs');

%% (all intermediates): sample pair from real data, go through all intermediates
% sample p(x3|x1) by going through each of the available xk. we sample through each xk a nSamples
% number of points. As before, we plot the detected covariance of samples p(x3|x1) vs the guessed
% covariance using our formula, but we also plot p(x3|x1) through each individual xk.

nTimes = 50;
nSamples = 1000;
nFeats = size(dataSigma, 1);

samplecov = zeros(nTimes, 1);
samplecovIndiv = zeros(nTimes, nFeats);
formulacov = zeros(nTimes, 1);
formulacov2 = zeros(nTimes, 1);
datacov = zeros(nTimes, 1);
for i = 1:nTimes

    % these are points 1,2,3 in our notation above. 2 is the intermediate we use.
    pair = randsample(nFeats, 2); 

    % sample
    [~, covSampleResult, loc, sigma] = sampleSigmaRecon(nSamples, dataSigma, pair);
    
    % estimate covariances
    samplecov(i) = covSampleResult(3); % sampling
    samplecovIndiv(i, :) = cellfun(@(x) x(3), loc)'; % sampling for each feature
    formulacov(i) = mean(sigma(:,2) .* sigma(:,5) ./ sigma(:,6)); % formula
    datacov(i) = dataSigma(pair(1), pair(2)); % real data
    
    % second formula (v2)
    sigma(pair, 2) = -inf;
    [~, mi] = max(sigma(:, 2)./sqrt(sigma(:, 6)));
    formulacov2(i) = sigma(mi, 2) .* sigma(mi, 5) ./ sigma(mi, 6);
end

% visualization
[~, si] = sort(samplecov);
figuresc(); hold on; ax = gca; ax.XMinorGrid = 'on';
trep = repmat((1:nTimes)', [1, nFeats]);
covSamplesIndivSorted = samplecovIndiv(si, :);
plot(trep(:), covSamplesIndivSorted(:), '.', sampleCovIndivColorSet{:});
plot(1:nTimes, samplecov(si), '*', sampleCovColorSet{:}); 
plot(1:nTimes, formulacov(si), 'o', formulaCovColorSet{:}); 
plot(1:nTimes, datacov(si), 'd', dataCovColorSet{:}); 
plot(1:nTimes, formulacov2(si), 'ok'); 
legend({'single-sample-$\sigma_{13}$', 'sample-$\sigma_{13}$', 'formula-$\sigma_{13}$', 'data-$\sigma_{13}$', 'formula-v2-$\sigma_{13}$'});
xlabel('runs (sorted by sample-$\sigma_{13}$)');
ylabel('estimated $\sigma_{13}$');
title('sigma estimates for different runs');

%% (all intermediates): use also the available weight covaraince of sparse real data
% instead of using true sigma, use iso-based model3 core sigma estimate (w/o reconstruction)
% this hsould help us see if the estimates/errors change by using the true first estimates. We
% shouldn't use the current reconstruction method since that's what we're trying to investigate
% here.
[~, m3dataSigmaIso] = model3sigma(isopatchesCen, maskpatches, 'none', 'none', nan);

% parameters
nTimes = 50;
nSamples = 1000;
nFeats = size(dataSigma, 1);

% sample the weights (maskSigma) from lowest to highest. 
% si will hold the indices into maskSigma for the trials
[~, si] = min(abs(bsxfun(@minus, maskSigma(:), linspace(1, max(maskSigma(:)), nTimes))), [], 1);

samplecov = zeros(nTimes, 1);
samplecovIndiv = zeros(nTimes, nFeats);
formulacov = zeros(nTimes, 1);
datacov = zeros(nTimes, 1);
for i = 1:nTimes

    % these are points 1, 2, 3 in our notation above. 2 is the intermediate we use.
    pair = ind2subvec([nFeats, nFeats], si(i));

    % build sigma and sample
    [~, covSampleResult, loc, sigma] = sampleSigmaRecon(nSamples, m3dataSigmaIso, pair);
    
    % estimate sigma
    samplecov(i) = covSampleResult(3); % sampling
    samplecovIndiv(i, :) = cellfun(@(x) x(3), loc)'; % individual feature sampling
    formulacov(i) = mean(sigma(:,2) .* sigma(:,5) ./ sigma(:,6)); % formula
    datacov(i) = dataSigma(pair(1), pair(2)); % estimate from real ISO data
    
    % using the second version of the formula. note, mi can be the current index i, if w_i is good.
    [~, mi] = max(sigma(:, 2) .* maskSigma(:, pair(1)) .* maskSigma(:, pair(2)) ./ sqrt(sigma(:, 6)));
    formulacov2(i) = sigma(mi, 2) .* sigma(mi, 5) ./ sigma(mi, 6);
end

% visualization 
wts = maskSigma(si);
wtsrep = repmat(wts(:), [1, nFeats]);
figuresc(); hold on; ax = gca; ax.XMinorGrid = 'on';
plot(wtsrep(:), samplecovIndiv(:), '.', sampleCovIndivColorSet{:});
plot(wts, samplecov, '*', sampleCovColorSet{:}); 
plot(wts, formulacov, 'o', formulaCovColorSet{:}); 
plot(wts, datacov, 'd', dataCovColorSet{:}); 
plot(wts, formulacov2, 'ok'); 
legend({'single-sample-$\sigma_{13}$', 'sample-$\sigma_{13}$ (via m3sigma)', 'formula-$\sigma_{13}$', 'data-$\sigma_{13}$', 'formula-v2-$\sigma_{13}$'});
xlabel('weight');
ylabel('estimated $\sigma_{13}$');
title('sigma estimates for different runs');

%% same as above, but sample from bisow model3 dataSigma
[~, m3dataSigmaBisow] = model3sigma(bisowpatchesCen, maskpatches, 'none', 'none', 0); 

% parameters
nTimes = 50;
nSamples = 1000;
nFeats = size(dataSigma, 1);

% sample the weights (maskSigma) from lowest to highest. 
% si will hold the indices into maskSigma for the trials
[~, si] = min(abs(bsxfun(@minus, maskSigma(:), linspace(1, max(maskSigma(:)), nTimes))), [], 1);

% sample
samplecov = zeros(nTimes, 1);
samplecovIndiv = zeros(nTimes, nFeats);
formulacov = zeros(nTimes, 1);
datacov = zeros(nTimes, 1);
for i = 1:nTimes

    % these are points 1, 2, 3 in our notation above. 2 is the intermediate we use.
    pair = ind2subvec([nFeats, nFeats], si(i));

    % sample
    [~, covSampleResult, loc, sigma] = sampleSigmaRecon(nSamples, m3dataSigmaBisow, pair); 
    
    % estimate covariances
    samplecov(i) = covSampleResult(3); % sampling
    samplecovIndiv(i, :) = cellfun(@(x) x(3), loc)'; % sampling individual features
    formulacov(i) = mean(sigma(:,2) .* sigma(:,5) ./ sigma(:,6)); % formula
    datacov(i) = dataSigma(pair(1), pair(2)); % real iso data
    
    % second formula
    [~, mi1] = max(sigma(:, 2) .* maskSigma(:, pair(1)) .* maskSigma(:, pair(2)) ./ sqrt(sigma(:, 6)));
    term1 = sigma(mi1, 2) .* sigma(mi1, 5) ./ sigma(mi1, 6);
    %[~, mi2] = max(sigma(:, 5) .* maskSigma(:, pair(1)) .* maskSigma(:, pair(2)) ./ sqrt(sigma(:, 6)));
    %term2 = sigma(mi2, 2) .* sigma(mi2, 5) ./ sigma(mi2, 6);
    formulacov2(i) = term1;
    %formulacov3(i) = (term1+term2)/2; % doesn't seem to help
    %formulacov4(i) = sqrt(term1.*term2); % doesn't seem to help
end

% visualization 
wts = maskSigma(si);
wtsrep = repmat(wts(:), [1, nFeats]);
figuresc(); hold on; ax = gca; ax.XMinorGrid = 'on';
plot(wtsrep(:), samplecovIndiv(:), '.', sampleCovIndivColorSet{:});
plot(wts, samplecov, '*', sampleCovColorSet{:}); 
plot(wts, formulacov, 'o', formulaCovColorSet{:}); 
plot(wts, datacov, 'd', dataCovColorSet{:}); 
plot(wts, formulacov2, 'ok'); 
legend({'single-sample-$\sigma_{13}$', 'sample-$\sigma_{13}$ (via m3sigma)', 'formula-$\sigma_{13}$', 'data-$\sigma_{13}$', 'formula-v2-$\sigma_{13}$', 'formula-v3-$\sigma_{13}$'});
xlabel('weight');
ylabel('estimated $\sigma_{13}$');
title('sigma estimates for different runs');

% TODO: also compute SSD for datacov vs formulacov, etc.

%% Some Lessons learned and things to try:
%  - somehow incorportate the interaction terms between xks?
%  whatever the first x2 doesn't explain of the covariance the second should
%  - a theoretical analysis of formula2. should you somehow use sigma(mi, 3)? etc.
