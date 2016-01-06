%% if cluster run
% TODO: should have common config files.

% break down to subvols.
% sgePrepareSubvols.sh trainmdfile, mods, volNames, patchSize, gridSpacing, atlVolSize, savefile

% run (via sge) training at each subvol
% sgeTrain.sh ...

% reconstruct gmm at each subvol

%% separate: cluster evaluation of K. This takes a long time, perhaps have to do very sparse

%% if solo run:
% prepare patch loading data
[loaddata, need] = papago.preloadTraining(gmmMethods, mods, opts.load, croprange);

% load subject data
[atlVolSize, subjmasks, dsregmaskvols, locVolumeAtlas, locVolumeSubjects] =  ...
    papago.preloadTesting(reconSubjs, mods, loadopts);


%% reconstruct patches in volume
[subgridlocs, subGridLength, subGridSize] = papago.prepgrid(atlVolSize, patchSize, patchOverlap, reconVolRange);
[fullgridlocs, fullGridLength, fullGridSize] = papago.prepgrid(atlVolSize, patchSize, 'sliding', reconVolRange);

% see vols2subvols which is currently being developed. 

% initiate reconstruction variables
reconSubs = cell(fullGridLength, 1);
reconPatches = cell(fullGridLength, 1);
weightrecon = cell(fullGridLength, 1);

% go through each point % TODO use verboseiter with text option.
for i = 1:subGridLength  % parfor
    loctic = tic;
    patches = [];
    
    % get gmm (load or learn)
    
   
    % initiate reconstruction
    fulloc = find(full2subloc == i);
    for f = fulloc(:)'
        reconSubs{f} = cell(numel(gmmMethods), 1);
        reconPatches{f} = cell(numel(gmmMethods), 1);
        weightrecon{f} = cell(numel(gmmMethods), 1);
    end
    
    % build the necessary gmm or wgmm
    gmms = cell(1, numel(gmmMethods));
    for t = 1:numel(gmmMethods)
        
        % train gmm
        gmms{t} = papago.train(gmmMethods{t}, patches, K, gmmargs, wgmmargs, reconSubj, pgmm{:});
        
        % impute using rotation
        % extract the corresponding patch from the subject volume provided
        volRange = arrayfunc(@(x, p) (x: x + p - 1)', location, patchSize);
        locPatchAtlas = cellfunc(@(ns) ns(volRange{:}), locVolumeAtlas);
        atlasPatch = dsregvol(volRange{:}); % same as dspathes(volidx == reconSubj & layeridx == 0, :)
        weightPatch = max(dsregmaskvol(volRange{:}), minW);
        
        m = 0;
        m = mean(atlasPatch(:)); % subtract mean;
        
        % try different reconstruction methods.
        % this needs a cleaner implementation of the conditional method, without for exampel the entire subject volume as input
        for f = fulloc(:)'
            [reconPatches{f}{t}, reconSubs{f}{t}] = ...
                papago.recon(gmms{t}, atlasPatch(:)' - m, exp(weightPatch(:)'-1), 'cond', ...
                'subjdsvol', subjdsvol - m, 'location', location, 'locPatchAtlas', locPatchAtlas, ...
                'locVolumeSubject', locVolumeSubject, 'subjmask', subjmask);
            reconPatches{f}{t} = reconPatches{f}{t} + m;
            weightrecon{f}{t} = papago.recon(gmms{t}, atlasPatch(:)' - m, exp(weightPatch(:)'-1), 'weig', 'percentile', 95) + m;
        end
        
        % a = atlasPatch; q = weigrecon{t}{i}; iii = isopatches(volidx == reconSubj & layeridx == 0, :); m = maskpatches(volidx == reconSubj & layeridx == 0, :);
        % view3Dopt(reshape(a, patchSize), reshape(q, patchSize), reshape(iii, patchSize), reshape(m, patchSize));
    end
    
    tc = toc(loctic);
    fprintf('%% done: %f %3.2fsec (total pred: %3.2fs)\n', i./subGridLength, tc, tc * subGridLength);
end

reconSubs1 = switchCellNesting(reconSubs);
reconPatches1 = switchCellNesting(reconPatches);
weightrecon1 = switchCellNesting(weightrecon);
