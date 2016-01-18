%% Test runtimes of gmms
% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
patchColPad = ones(1, 3) * 2;

regVal = 1e-4; % regulaization to the diagonal of subjSigma, if using method forward

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load *ground truth* data column from buckner
% bucknermd = loadmd([SYNTHESIS_DATA_PATH, 'buckner', '_restor_md_*']);
% [bucknerIsoPatchCol] = ...
%    subspacetools.md2patchcol(bucknermd, 'brainIso2Ds5Us5size', atlPatchSize, atlLoc, patchColPad);

%% compute gaussian mixture model
% learn a gaussian mixture model with K clusters from the *true* isotropic data, 

% compute the gaussian mixture model
nRepls = 3;
tolfuns = [0.1, 0.01, 0.001, 0.0001];
gmmks = [1, 2, 3, 5, 10, 15, 25];
times = zeros(numel(gmmks), numel(tolfuns));
for ki = 1:numel(gmmks)
    k = gmmks(ki);
    for ti = 1:numel(tolfuns)
        t = tolfuns(ti);
        gmmopt = statset('Display', 'off', 'MaxIter', 100, 'TolFun', t);
        
        tic
        gmm = fitgmdist(bucknerIsoPatchCol, k, 'regularizationValue', regVal, ...
            'replicates', nRepls, 'Options', gmmopt);
        times(ki, ti) = toc;
        fprintf(2, 'GMM (%d replicate(s), k=%d, tol=%f), took %3.3f sec\n', ...
            nRepls, k, t, times(ki, ti));
    end
end

%% plot times.
plot(gmmks, times, '.-');
names = arrayfunc(@(x) sprintf('tol %5.4f', x), tolfuns);
legend(names);
xlabel('nr clusters');
ylabel('time (sec)');
grid on;

