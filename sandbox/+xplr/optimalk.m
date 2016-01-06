N = 5;
krange = [1, 2, 3, 5, 10, 25];

% TODO: need to load data as in testVolIsoGMMRecons_Core.

% Setup sparse (ish)
[subgridlocs, subGridLength, subGridSize] = subspacetools.prepgrid(atlVolSize, patchSize, 'discrete');

for i = 1:subGridLength  % parfor
    loctic = tic;
    patches = [];

    % get location of this patch
    location = subgridlocs(i, :);
    
    % load data
    [p, patches.layeridx, patches.volidx] = subspacetools.loadpatches(loadmethod, loaddata{:}, ...
        patchSize, location, locpad);
    if needisovols || needbisowvols, patches.iso = p{1};  p = p(2:end); end
    if needmasks, patches.mask = p{1}; patches.W = max(patches.mask, minW); p = p(2:end); end
    if needdsvols, patches.ds = p{1}; p = p(2:end); end
    if needbisowvols, patches.ds = patches.mask .* patches.iso + (1 - patches.mask) .* p{2}; end
    


    %krange = linspace(1, 8, 8);


    patches.iso = bsxfun(@minus, patches.iso, mean(patches.iso, 2));
    trainidx = false(size(patches.iso, 1), N);
    for n = 1:N
        trainind = randsample(size(patches.iso, 1), round(size(patches.iso, 1)/2));
        trainidx(trainind, n) = true;
    end



    % prepare variables.
    trainlogp = zeros(numel(krange), N);
    testlogp = zeros(numel(krange), N);
    timet = zeros(numel(krange), N);


    parfor ki = 1:numel(krange)
        k = krange(ki);
        trainp = [];
        testp = []

        for n = 1:N
            trainp.iso = patches.iso(trainidx(:, n), :);
            testp.iso = patches.iso(~trainidx(:, n), :);

            tic;
            isogmm = papago.train('iso-gmm', trainp, k, gmmargs, wgmmargs, reconSubj, pgmm{:});
            isogmm.logdetw = zeros(size(trainp.iso, 1), 1); % wgmm.logdet(diag(W(i, :)))
            timet(ki, n) = toc;

            trainlogp(ki, n) = isogmm.logp(trainp.iso, trainp.iso*0+1);
            testlogp(ki, n) = isogmm.logp(testp.iso, testp.iso*0+1);

            
            
            [reconPatches{i}{t}, reconSubs{i}{t}] = ...
                papago.recon(gmms{t}, atlasPatch(:)' - m, exp(weightPatch(:)'-1), 'cond', ...
                'subjdsvol', subjdsvol - m, 'location', location, 'locPatchAtlas', locPatchAtlas, ...
                'locVolumeSubject', locVolumeSubject, 'subjmask', subjmask);
            
        end
        ki
    end

    %% visualize
    figure(); 
    ttl = sprintf('location: %d %d %d', location);
    subplot(1, 2, 1);
    errorbar(krange, mean(trainlogp, 2), std(trainlogp, [], 2)); hold on; axislabels('K', 'loglik', ttl);
    errorbar(krange, mean(testlogp, 2), std(testlogp, [], 2)); legend({'train', 'test'});
    subplot(1, 2, 2);
    errorbar(krange, mean(timet, 2), std(timet, [], 2)); axislabels('K', 'time taken', ttl);
    
    drawnow;
    
    [~, mi] = max(mean(testlogp, 2));
    optK(i) = krange(mi);
end



z = nan(atlVolSize);
z(subvec2ind(atlVolSize, subgridlocs)) = optK;
zi = inpaintn(z);
view3Dopt(nii2vol(loadNii(BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS5_US5)), z, zi);

