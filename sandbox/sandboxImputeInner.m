%% TODO: combine with sandboxImpute

%% setup
% md should be loaded from patchstartup.m
patchSize = ones(1, 3) * 7;
location = [56, 72, 92];
locpad = ones(1, 3) * 2;
locsearch = [7, 7, 1]; %ones(1, 3) * 3;

% load datasets
niis = md2niis(md, dataset);

%% Trial 1
location = [56, 72, 92];
locsearch = [3, 3, 0];
prange = arrayfunc(@(x, d, p) (x - d: x + d + p - 1)', location, locsearch, patchSize);
resop = @(x) volresize(x.img, round(size(x.img)));

vols = cellfunc(@(x) x(prange{:}), cellfunc(resop, niis.ds));
validvols = cellfunc(@(x) x(prange{:}) > 0.1, cellfunc(resop, niis.mask));
isovols = cellfunc(@(x) x(prange{:}), cellfunc(resop, niis.iso));


%%
quiltedvolumes = volPatchImpute(vols, validvols, patchSize, [1, 1, 1], isovols);

%%
[recovpatches, meanX] = pcaimpute(tmppatches, K, nRep, isopatches);