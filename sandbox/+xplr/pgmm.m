%% Explore global GMM (GPGMM and GGMM) 
import subspacetools.reshapeN3Dto2D;

% load pgmm and ggmm
load(fullfile(SYNTHESIS_DATA_PATH, 'online_gmm_million_patches__'), 'gmdistfull', 'gmdistviapca');
ggmm = gmdistfull;
gpgmm = gmdistviapca;

%% location-specific.
location = [31, 31, 21];
[patches.iso, patches.layeridx, patches.volidx] = subspacetools.loadpatches(loadmethod, mfptrs, isoregmod, ...
    patchSize, location, locpad);
patches.iso = bsxfun(@minus, patches.iso, mean(patches.iso, 2));

% load gmms at location. high r/i, or low r/i
pgmm = {};
tic
gmmargs = {'RegularizationValue', regval, 'replicates', 3, 'Options', statset('Display', 'final', 'MaxIter', 5, 'TolFun', 0.01)};
gmmlocr3i5 = papago.train('iso-gmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
gmmlocr3i5toc = toc(); tic
gmmargs = {'RegularizationValue', regval, 'replicates', 1, 'Options', statset('Display', 'final', 'MaxIter', 1, 'TolFun', 0.01)};
gmmlocr1i1 = papago.train('iso-gmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
gmmlocr1i1toc = toc(); tic

% load gpgmm and ggmm at location. direct locat adaptation, or plus low r/i
pgmm = {gpgmm};
gmmargs = {'RegularizationValue', regval, 'replicates', 3, 'Options', statset('Display', 'final', 'MaxIter', 5, 'TolFun', 0.01)};
gpgmmloc = papago.train('iso-pgmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
gpgmmloctoc = toc(); tic;
S = gmm2Start(gpgmmloc, 'wparameters');
gmmargs = {'RegularizationValue', regval, 'replicates', 1, 'Options', statset('Display', 'final', 'MaxIter', 1, 'TolFun', 0.01), 'Start', S};
gpgmmlocr1i1 = papago.train('iso-gmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
gpgmmlocr1i1toc = toc(); tic

pgmm = {ggmm};
gmmargs = {'RegularizationValue', regval, 'replicates', 3, 'Options', statset('Display', 'final', 'MaxIter', 5, 'TolFun', 0.01)};
ggmmloc = papago.train('iso-pgmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
ggmmloctoc = toc(); tic
S = gmm2Start(ggmmloc, 'wparameters');
gmmargs = {'RegularizationValue', regval, 'replicates', 1, 'Options', statset('Display', 'final', 'MaxIter', 1, 'TolFun', 0.01), 'Start', S};
ggmmlocr1i1 = papago.train('iso-gmm', patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
ggmmlocr1i1toc = toc(); tic


gmms = {gmmlocr3i5, gmmlocr1i1, gpgmmloc, gpgmmlocr1i1, ggmmloc, ggmmlocr1i1};
tocs = {gmmlocr3i5toc, gmmlocr1i1toc, gpgmmloctoc, gpgmmlocr1i1toc, ggmmloctoc, ggmmlocr1i1toc};
gmmTitles = {'gmmlocr3i5', 'gmmlocr1i1', 'gpgmmloc', 'gpgmmlocr1i1', 'ggmmloc', 'ggmmlocr1i1'};
gmmTitles = cellfunc(@(g, t) sprintf('%s %1.3f sec', g, t), gmmTitles, tocs);
compareGMMmeans(gmms, patchSize, gmmTitles);

%% Location GMM vs PGMM

% assumes we have a gmm computed manually at a location (gmmloc), a general pgmm (pgmm) , the pgmm
% re-evaluated at that location (pgmmloc);


compareGMMs({gmmloc, pgmm, pgmmloc}, patchSize, {'gmm@loc mus', 'pgmm mus', 'pgmm@loc mus'});