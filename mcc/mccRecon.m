function mccRecon(gmmfile, testmdfile, modalities, subjidx, subtractpatchmean, minW, patchSize, subjoutfile)
% for one subject, reconstruct all patches sub this subvolumes using the given gmm.
%
% take in (md, subjidx, matmodality, modalities, gmmfile, subtractpatchmean, patchSize, subjoutfile)
%
%
%
% steps:
% 1. load subject files


    % modalities needs to be in this order: ds, dsmask, dsreg, dsregmask
    if ischar(modalities), modalities = str2cell([modalities ','], ','); end
    mods.ds = modalities{1};
    mods.dsmask = modalities{2};
    mods.dsreg = modalities{3};
    mods.dsregmask = modalities{4};
    
    % prep data
    [atlVolSize, subjmasks, dsregmaskvols, locVolumeAtlas, locVolumeSubjects, dsregvols] =  ...
        preloadTesting(subjidx, mods, testmdfile);

    % get data location
    q = load(subvolfile, 'nfo');
    nfo = q.nfo;
        
    % load gmm
    load(gmmfile);
    
    % need:
    dsregvol = dsregvols{1};
    dsregmaskvol = dsregmaskvols{1};
    
    % setup locations
    subs = patchlib.grid(nfo.subvolSize, patchSize, 'sub');
    
    reconPatches = zeros(size(subs, 1), prod(patchSize));
    reconSubs = zeros(size(subs, 1), numel(patchSize));
    for i = 1:size(subs, 1)
        location = nfo.location + subs(i, :);

        % impute using rotation
        % extract the corresponding patch from the subject volume provided
        volRange = arrayfunc(@(x, p) (x: x + p - 1)', location, patchSize);
        locPatchAtlas = cellfunc(@(ns) ns(volRange{:}), locVolumeAtlas);
        atlasPatch = dsregvol(volRange{:}); 
        weightPatch = max(dsregmaskvol(volRange{:}), minW);

        % subtract means if necessary
        if subtractpatchmean
            patchMean = mean(atlasPatch);
            atlasPatch = atlasPatch(:)' - patchMean;
        end

        [reconPatches(i, :), reconSubs(i, :)] = papago.recon(gmm, atlasPatch, weightPatch, 'cond', ...
            'subjdsvol', subjdsvol - patchMean, 'location', location, 'locPatchAtlas', locPatchAtlas, ...
            'locVolumeSubject', locVolumeSubject, 'subjmask', subjmask);
        
        if subtractpatchmean
            reconPatches(i, :) = reconPatches(i, :) + patchMean;
        end

    end
        
    % save reconPatches.
    save(subjoutfile, 'reconPatches', 'reconSubs');
    