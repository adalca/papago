function [loaddata, need] = preloadTraining(gmmMethods, mods, loadopts, croprange)
% PRELOADTRAINING
%   [loaddata, need] = preloadTraining(gmmMethods, mods, loadopts) preload training data (volumes or
%   matfiles) for papago. loaddata is a cell array that can be passed into loadpatches.
%
%   [loaddata, need] = preloadTraining(gmmMethods, mods, loadopts, croprange) allows the
%   specification of a cropping range. croprange is a 1xnDims cell, each entry specifying the rows
%   to be kept. e.g. croprange = {1:10, 2:50, 10:30}. only used with the 'volumes' load method.
%
%   gmmMethods - a cell array of gmm methods for which we need to load training data.
%   mods - a struct with the name of the following modalities: 
%       isoreg, dsmaskreg, dsreg, bisoreg
%   loadopts - a struct with parameters:
%       method - 'volumes' or 'matfiles'
%       atldstype - 'ds' or 'bisow' --- the type of "downsampled" data
%       trainmdfile, 
%       trainmfptrsfile
%   <croprange> - optional - the range to crop volumes to.

    % prepare loading data
    load(loadopts.trainmfptrsfile); 
    load(loadopts.trainmdfile); 

    % some logicals determining what we need
    need.isovols = any(strncmp('iso', gmmMethods, 3));
    dsdata = any(strncmp('ds', gmmMethods, 2));
    wgmmmethod = any(cellfun(@(x) ~isempty(x), strfind(gmmMethods, 'wgmm')));
    need.masks = dsdata || wgmmmethod;
    need.dsvols = dsdata && strcmp(loadopts.atldstype, 'ds');
    need.bisowvols = dsdata && strcmp(loadopts.atldstype, 'bisow');

    % prepare volumes if necessary
    tic
    if strcmp(loadopts.method, 'volumes'); 
        % use md instead of mfptrs for volumes. It's much faster. 
        % TODO: Why? mfptrs loading not well implemented?
        if exist('croprange', 'var')
            getpfn = @(mod) subspacetools.loadpatches('md', md, mod, @(x) x(croprange{:})); 
        else
            getpfn = @(mod) subspacetools.loadpatches('md', md, mod); 
        end

        loaddata = {{}};
        if need.isovols || need.bisowvols; loaddata{1} = [loaddata{1}, getpfn(mods.isoreg)]; end
        if need.masks; loaddata{1} = [loaddata{1}, getpfn(mods.dsmaskreg)]; end
        if need.dsvols; loaddata{1} = [loaddata{1}, getpfn(mods.dsreg)]; end
        if need.bisowvols; loaddata{1} = [loaddata{1}, getpfn(mods.bisoreg)]; end
        % bisowvols = cellfunc(@(m, i, b)  m .* i + (1-m) .* b, maskvols, isovols, bisovols);

    % prepare matfiles if necessary
    else % 
        assert(strcmp(loadopts.method, 'matfiles'));
        
        loaddata = {mfptrs, {}};
        if need.isovols || need.bisowvols, loaddata{2} = [loaddata{2}, mods.isoreg]; end
        if need.masks, loaddata{2} = [loaddata{2}, mods.dsmaskreg]; end
        if need.dsvols, loaddata{2} = [loaddata{2}, mods.dsreg]; end
        if need.bisowvols, loaddata{2} = [loaddata{2}, mods.bisoreg]; end
    end
    fprintf('Prepared volumes in %1.3f\n', toc);
    