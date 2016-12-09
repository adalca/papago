function subvolReconFlexible(wgmmFile, atlSubvols, iniReconFile, subjVol, subjMask, interpMatData, subjoutFile)
%code for LSR recon:
%
% for each subvol
%   get subvol data using vol2subvolInterpData
%   break into R pieces using code from loadData (?)
%   call wgmm.recon.
%   output the subvol.

    % load in data
    if ischar(subjVol) && sys.isfile(subjVol), subjVol = nii2vol(subjVol); end
    if ischar(subjMask) && sys.isfile(subjMask), subjMask = nii2vol(subjMask); end
    if ischar(interpMatData) && sys.isfile(interpMatData), interpMatData = load(interpMatData); end
     
    % atlSubvols could be a list. if it is, recursively call subvolReconFlexible. Note this has to
    % be after volumes were loaded, or there's no help.
    if ~sys.isfile(atlSubvols)
        % assume it's a list
        atlFiles = strsplit(atlSubvols, ',');
        outFiles = strsplit(atlSubvols, ',');
        wgmmFiles = strsplit(wgmmFile, ',');
        
        assert(numel(atlFiles) == numel(outFiles), ...
            'The number of subvolume files does not match the number of output files');
        assert(numel(atlFiles) == numel(wgmmFiles), ...
            'The number of subvolume files does not match the number of wgmm files');
        for fi = 1:numel(atlFiles)
            subvolReconFlexible(wgmmFiles{fi}, atlFiles{fi}, iniReconFile, subjVol, subjMask, interpMatData, outFiles{fi});
        end
    end
    
    % load subvolume location
    q = load(atlSubvols);
    reconLoc = q.nfo.subvolLoc;
    atlSubvolSize = q.nfo.subvolSize;
    
    % get other ingo.
    ini = ini2struct(iniReconFile);
    atlPatchSize = ini.atlPatchSize;
    reconModel = ini.reconModel;
    
    % load volumes
    q = load(wgmmFile);
    wg = q.wg;
    
    % get rotation data for this subvolume. do NOT pass in an output file -- we just want this
    % returned. This loads a signle volume, despite the plural (e.g. subvols) in the variable names.
    [interpData.subVols, interpData.subVolMasks, interpData.atl2SubjInterpMat, interpData.subj2AtlInterpMat] = ...
        vol2subvolInterpData(interpMatData, subjVol, subjMask, reconLoc, atlSubvolSize);
    
    % get R matrices.
    interpDataPatches = interpData2patchRotData(interpData, subvolSize, 0, atlPatchSize, 1);

    % reconstruct patches
    nPatches = numel(interpDataPatches.yorigAll);
    volRotData = paffine.prepareWgmmRData(interpDataPatches, ones(1, nPatches));
    reconPatches = wg.recon(volRotData, reconModel); %'latentMissingR');
    
    % reconstruct subvolume
    rotmasks = interpDataPatches.yrotmasksAll(:);
    patches = cellfunc(@maskvox2vol, reconPatches(:), rotmasks);
    recmasks = cellfunc(@maskvox2vol, volRotData.rWeight, rotmasks);
    [reconVol, reconWeight] = patchlib.quiltIrregularPatches(volStarts, patches, ...
        'volSize', atlSubvolSize, 'weightPatches', recmasks);

    % save reconPatches.
    save(subjoutFile, 'reconVol', 'reconLoc', 'reconWeight');
