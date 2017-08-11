function mccSubvolReconFlexibleWrap(wgmmFileFormat, atlSubvolsFormat, iniReconFile, ...
    atlSubjVol, atlSubjMask, interpMatData, subjoutFileFormat, subvolIdsFile)

    % load stuff you only need to load in once - R file, subject files
    tic; fprintf('Reading in common batch files...');
    if ischar(atlSubjVol) && sys.isfile(atlSubjVol), atlSubjVol = nii2vol(atlSubjVol); end
    if ischar(atlSubjMask) && sys.isfile(atlSubjMask), atlSubjMask = nii2vol(atlSubjMask); end
    if ischar(interpMatData) && sys.isfile(interpMatData), interpMatData = load(interpMatData, 'atlLoc2SubjSpace', 'atl2subjR'); end
    fprintf(' done (%3.2f)\n', toc);
    
    % load in subvol ids
    tic; fprintf('Reading in subvolIds...');
    if sys.isfile(subvolIdsFile)
        fid = fopen(subvolIdsFile);
        subvolIds = textscan(fid, '%d\n');
        fclose(fid);
        subvolIds = cat(1, subvolIds{:});
    else
        error('subvolIdsFile for now has to be a file');
    end
    fprintf(' done. Found %d ids in (%3.2f)\n', numel(subvolIds), toc);
    
    for vi = 1:numel(subvolIds)
        subvolId = subvolIds(vi);
        tic; fprintf('\nStarting subvol %d\n', subvolId);
        
        % get the wgmm, subvol and output files
        wgmmFile = sprintf(wgmmFileFormat, subvolId);
        atlSubvols = sprintf(atlSubvolsFormat, subvolId);
        subjoutFile = sprintf(subjoutFileFormat, subvolId);
        
        fprintf('wgmmFile: %s\n', wgmmFile);
        fprintf('atlSubvols: %s\n', atlSubvols);
        fprintf('subjoutFile: %s\n\n', subjoutFile);
        
        
        
        % call recon
        try
            subvol.subvolReconFlexible(wgmmFile, atlSubvols, iniReconFile, ...
                atlSubjVol, atlSubjMask, interpMatData, subjoutFile);
        catch err
            err;
        end
        
        fprintf('Done subvolid %d in %3.2f\n', subvolId, toc);
        dispDashedLine;
    end
end