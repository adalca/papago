%% explore an online gmm
error('unfinished');

% parameters
patchSize = ones(1, 3) * 9;
K = 50; 
maxOGMMIters = 100;
gmmopt = statset('Display', 'iter', 'MaxIter', 25);
gmmregval = 1e-4;
weightThresh = 1e-7;
 
% online parameters
nInitialPatches = 1000000;
initmethod = 'pca-wgmm'; % 'standard-wgmm', 'pca-wgmm', % also inits with gmm?
startmethod = 'assignment'; % 'assignment' or 'parameters'
initvarExplained = 99; % for 'pca-gmm'
nStepSamples = 20000; 

% loading and saving parameters
% loadargs = {subspacetools.loadpatches('niistruct', 'data/bucknerNiis_5_5.mat', 'niis', 'iso')};
loadargs = {subspacetools.loadpatches('md', 'data/adni_md.mat', 'matfile', 'rigidRegBrain'), 'rigidRegBrain', 'regBrainDs5Us5Mask'};
savepath = [DATA_PATH, 'output/subspace/onlinegmm/'];

%% Online GMM initialization

% get initial data
[patches, locsamples, volsamples] = patchlib.vol2samples(nInitialPatches, patchSize, loadargs{:});

% initialize gmm
switch initmethod
    case 'standard-wgmm'
        gmdist = fitgmdist(patches, K, 'Replicates', 3, 'RegularizationValue', gmmregval, 'Options', gmmopt);
        gminit = gmdist;
        S = subspacetools.gmm2Start(gmdist, 'parameters');
        
    case 'pca-wgmm'
        [coeff, score, ~, ~, explained, mu] = pca(patches);
        cexplained = cumsum(explained);
        f = find(cexplained > initvarExplained, 1, 'first');
        fprintf('%d components explain %3.1f%% of variance\n', f, cexplained(f));
        gmdist = fitgmdist(score(:, 1:f), K, 'Replicates', 3, 'RegularizationValue', gmmregval, 'Options', gmmopt);
        
        gminit = gmdist;
        if strcmp(starttype, 'parameters')
            % get parameters by doing a single iteration of a full-gmm
            gmdist = fitgmdist(patches, K, 'Replicates', 1, 'Regularize', gmmregval, ...
                'Options', statset('Display', 'iter', 'MaxIter', 1), 'Start', S);
        end
        
    otherwise
        error('onlineGMM: Unknown init method');
end

robustSave(sprintf('%scurrentGMM_iter0.mat', savepath), 'gmdist', 'gminit', 'locsamples', 'volsamples');

%% Online GMM updates
gmmopt = statset('Display', 'iter', 'MaxIter', 7);
gmmpar = {'Replicates', 1, 'RegularizationValue', gmmregval, 'Options', gmmopt};

% update steps
for i = 1:maxOGMMIters

    % get more data!
    [patches, locsamples, volsamples] = patchlib.vol2samples(nStepSamples, patchSize, loadargs{:});
    
    % prepare 'Start' for next GMM. 
    % If it's the first iteration and 'assignment' prepare pca-projected data.
    data = patches;
    if i == 1 && strcmp(startmethod, 'assignment')
        data = pcax.project(bsxfun(@minus, patches, mu)', coeff(:, 1:f))';
    end
    S = subspacetools.gmm2Start(gmdist, startmethod, data);
    
    % update
    gmdist = fitgmdist(patches, K, gmmpar{:}, 'Start', S);
    [~, nloglik] = gmdist.posterior(patches);
    
    % save
    if (mod(i, 10) == 0)
        robustSave(sprintf('%scurrentGMM_iter%d.mat', savepath, i), ...
            'gmdist', 'nloglik', 'locsamples', 'volsamples');
    end
    
    fprintf('oiter: %d. nloglik:%f\n', i, nloglik); 
end
