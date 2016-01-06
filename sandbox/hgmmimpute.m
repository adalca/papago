function [vols, Xiw, gmdw] = hgmmimpute(niifiles, dsRate, usRates, patchFactors, location, locpad, nClust, isoniis, K, nReps, reconmethod)
    o3 = ones(1, 3);

    for u = 1:numel(usRates)
        usRate = usRates(u);

        % compute ideal location
        sclocation = round(location * usRate./dsRate); 
        scPatchSize = o3 * patchFactors(u);
        scerrfun = @(x, y) subspacetools.patcherror(x, y, @msd, 'gauss', scPatchSize, 1.5);

        % load appropriate niftis
%             scniis = subspacetools.md2niis(md, 'buckner', dsRate, usRate);    
        scniisx = load(niifiles{u});
        scniis = scniisx.niis;
        % extract appropriate patches at this level
        [dspatches, dspatchidx] = subspacetools.nii2patchcol(scniis.ds, scPatchSize, sclocation, locpad);
        [maskpatches, maskpatchidx] = subspacetools.nii2patchcol(scniis.mask, scPatchSize, sclocation, locpad);

        % extract previous volume's patch at this location
        if u == 1
            prevvols = cellfunc(@(x) x.img, scniis.ds);
        else
            % TODO: check how upsampling and downsam,ping is done so we align.
            % Do same method of upsampling or make sure it fits!?
            prevvols = cellfunc(@(x, y) volresize(x, size(y.img)), vols, scniis.ds);
        end
        [prevpatches, prevpatchidx] = subspacetools.vols2patchcol(prevvols, scPatchSize, sclocation, locpad);

        % combine patches via masks
        patches = dspatches .* maskpatches + prevpatches .* (1 - maskpatches);

        % perform gmm reconstruction (using weights)
        switch reconmethod
            case 'method1'
                fakeTrueVols = cellfunc(@(x, y) volresize(x, size(y)), isoniis, scniis.ds);
                [fisopatches, fisopatchidx] = subspacetools.nii2patchcol(fakeTrueVols, scPatchSize, sclocation, locpad);
                fprintf('Linear Interp: %e\n', scerrfun(dspatches, fisopatches));
                opts = {'valinit', patches, 'trueData', fisopatches, 'realErrFun', scerrfun};
                [Xiw, gmdw] = gmmimpute(patches, nClust, K, 'maxIter', nReps, opts{:}, 'valwts', maskpatches);
            case 'model3-rotrecon'
                
           
            otherwise
                error('Unknown GMM learn and patch reconstruction method');
        end
            

        % patches to volume.
        patchsel = Xiw(dspatchidx, :);
        vols = prevvols;
        r = arrayfunc(@(x, p) x:x+p-1, sclocation, scPatchSize);
        for i = 1:numel(vols)
            vols{i}(r{:}) = reshape(patchsel(i, :), scPatchSize);
        end
    end
end