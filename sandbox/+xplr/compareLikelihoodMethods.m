%%  Compare Likelihood (and psoterior) patch assignment methods given GMM for Sparse data
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 
% minus this subject.

%% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
gmmK = 5; % 5 good for ventricle, 25 for cortex?
crmethod = 'inverse'; % 'forward', 'inverse'
regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward
reconSubj = 3; %1, 3
patchColPad = ones(1, 3) * 2;

% settings
nCompareReps = 100;
nLogMethods = 4;

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load *ground truth* data column from buckner
bucknermd = loadmd([SYNTHESIS_DATA_PATH, 'buckner', '_restor_md_*']);
[bucknerIsoPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(bucknermd, 'brainIso2Ds5Us5size', atlPatchSize, atlLoc, patchColPad);

% load selected ADNI subject volumes
adnimd = loadmd([SYNTHESIS_DATA_PATH, 'adni', '_restor_md_*']);
dsSubjNii = adnimd.loadModality('brainDs5Us5', reconSubj);
subjVol = double(dsSubjNii.img);
subjWeightVol = logical(adnimd.loadVolume('brainDs5Us5Mask', reconSubj));
dsSubjInAtlNii = adnimd.loadModality('brainDs5Us5Reg', reconSubj);
subjInAtlTform = load(adnimd.getModality('brainDs5Us5RegMat', reconSubj));

% prepare necessary inputs for conditional-based reconstruction
subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
tform = subjInAtlTform.tform;
atlVolSize = size(dsSubjInAtlNii.img);
subjLoc2AtlSpace = tform2cor3d(tform, size(subjVol), subjDims, atlVolSize, atlDims);
atlLoc2SubjSpace = tform2cor3d(tform, size(subjVol), subjDims, atlVolSize, atlDims, 'backward');
extraReconArg = ifelse(strcmp(crmethod, 'inverse'), subjLoc2AtlSpace, regVal);
dsSubjInAtlMasknii = bucknermd.loadModality('brainDs5Us5RegMask', reconSubj);

%% patches
isoSubjVol = adnimd.loadVolume('brainCropped5', reconSubj);

% atlas patches
dsAtlPatch = cropVolume(dsSubjInAtlNii.img, atlLoc, atlLoc + atlPatchSize - 1);
dsAtlMaskPatch = cropVolume(dsSubjInAtlMasknii.img, atlLoc, atlLoc + atlPatchSize - 1);

% subject patches
[~, subjInterpMask, subjPatchLoc, subjPatchSize] = ...
    paffine.volCor2patchInterpmat(method, atlLoc, atlPatchSize, atlLoc2SubjSpace, subjLoc2AtlSpace);

dsSubjPatch = cropVolume(subjVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
correctSubjPatch = cropVolume(isoSubjVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
maskSubjPatch = cropVolume(subjWeightVol, subjPatchLoc, subjPatchLoc + subjPatchSize - 1);
atlInSubjMask = subjInterpMask; 

%% run GMM and logs

% gmm options
gmmopt = statset('Display', 'off', 'MaxIter', 20, 'TolFun', 0.001);

% prepare stats
stats.cntv = zeros(4, nCompareReps);
stats.rankv = zeros(4, nCompareReps);
stats.percv = zeros(4, nCompareReps);

% save the gmms
gmm = cell(1, nCompareReps);

fprintf('%7s | %10s | %15s | %20s | %10s | %15s | %20s\n', ..., 
    'itr', 'count correct', 'rank', 'ssd diff to best', ...
    'sum of counts', 'avg rank', 'mean sdd diff to best');
for rep = 1:nCompareReps
    % compute the gaussian mixture model
    gmm{rep} = fitgmdist(bucknerIsoPatchCol, gmmK, ...
        'regularizationValue', regVal, 'replicates', 5, 'Options', gmmopt);

    % reconstruct
    logp = repmat({zeros(1, K)}, [1, nLogMethods]);
    reconPatches = cell(1, K);
    subjPatchLocs = cell(1, K);
    
    for k = 1:K
        atlMu = gmm{rep}.mu(k, :)';
        atlSigma = gmm{rep}.Sigma(:, :, k);
       
        % v1: compute logp via subspace of covariance representing the known voxels in subject space
        % also compute the reconstruction using this method
        [reconPatches{k}, subjPatchLocs{k}, logp{1}(k)] = paffine.recon(atlMu, atlSigma + eye(prod(atlPatchSize)) * regVal, atlLoc, atlPatchSize, ...
            subjVol, subjWeightVol, atlLoc2SubjSpace, method, extraReconArg);
        
        % v2: compute logp via whole volume, weighted.
        [subjMu, subjSigma] = ...
            paffine.atl2SubjGauss(atlMu, atlSigma + eye(prod(atlPatchSize)) * regVal, method, atlLoc, atlPatchSize, atlLoc2SubjSpace, extraReconArg);
        x = dsSubjPatch(:) .* maskSubjPatch(:); x(~subjInterpMask) = [];
        mu = subjMu(:) .* maskSubjPatch(:); mu(~subjInterpMask) = [];
        sigma = subjSigma(subjInterpMask, subjInterpMask);
        logp{2}(k) = logmvnpdf(x(:)', mu(:)', sigma);
        
        % v3: don't include weights.
        x = dsSubjPatch(:); x(~subjInterpMask) = [];
        mu = subjMu(:); mu(~subjInterpMask) = [];
        sigma = subjSigma(subjInterpMask, subjInterpMask);
        logp{3}(k) = logmvnpdf(x(:)', mu(:)', sigma);
        
        % v4: atlas space of v2
        logp{4}(k) = logmvnpdf(dsAtlPatch(:)' .* dsAtlMaskPatch(:)', atlMu(:)' .* dsAtlMaskPatch(:)', atlSigma);
    end
    
    % compute true ssd distance
    ssds = cellfun(@(x) ssd(x(atlInSubjMask), correctSubjPatch(atlInSubjMask)), reconPatches);
    [~, besti] = min(ssds);
    
    % compute measures
    for i = 1:length(logp), 
        % compute min logs
        meas = log(gmm{rep}.ComponentProportion) + logp{i};
        [~, bestlogpv] = max(meas);
    
        % 0-1 measure
        stats.cntv(i, rep) = bestlogpv == besti;
        
        % ranking 
        [~, locrankv] = sort(logp{i}, 'descend');
        stats.rankv(i, rep) = find(locrankv == besti);
        
        % percent of logp
        % percv(rep) = 1 - exp(logp1(besti) - max(logp1));
        
        % percent of ssd
        stats.percv(i, rep) = (ssds(bestlogpv) - ssds(besti)) ./ ssds(besti);
    end
    
    % print this rep's stats
    fprintf('%5d\t', rep);
    fprintf(' | '); fprintf('%3d ', stats.cntv(:, rep));
    fprintf(' | '); fprintf('%5.1f ', stats.rankv(:, rep));
    fprintf(' | '); fprintf('%7.2f ', stats.percv(:, rep));
    
    % print summary stats up to and including this rep
    fprintf(' | '); fprintf('%3d ', sum(stats.cntv(:, 1:rep), 2));
    fprintf('| '); fprintf('%5.1f ', mean(stats.rankv(:, 1:rep), 2));
    fprintf('| '); fprintf('%7.2f ', mean(stats.percv(:, 1:rep), 2));
    newline;
end

