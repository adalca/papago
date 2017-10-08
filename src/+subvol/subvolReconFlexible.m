function subvolReconFlexible(wgmmFile, atlSubvols, iniReconFile, atlSubjVol, atlSubjMask, interpMatData, subjoutFile)
%code for LSR recon:
%
% for each subvol
%   get subvol data using vol2subvolInterpData
%   break into R pieces using code from loadData (?)
%   call wgmm.recon.
%   output the subvol.

    % load in data
    if ischar(atlSubjVol) && sys.isfile(atlSubjVol), atlSubjVol = nii2vol(atlSubjVol); end
    if ischar(atlSubjMask) && sys.isfile(atlSubjMask), atlSubjMask = nii2vol(atlSubjMask); end
    if ischar(interpMatData) && sys.isfile(interpMatData), interpMatData = load(interpMatData); end
     
    % atlSubvols could be a list. if it is, recursively call subvolReconFlexible. Note this has to
    % be after volumes were loaded, or there's no help.
    if ~sys.isfile(atlSubvols)
        % assume it's a list
        atlFiles = strsplit(atlSubvols, ',');
        outFiles = strsplit(subjoutFile, ',');
        wgmmFiles = strsplit(wgmmFile, ',');
        
        assert(numel(atlFiles) == numel(outFiles), ...
            'The number of subvolume files does not match the number of output files');
        assert(numel(atlFiles) == numel(wgmmFiles), ...
            'The number of subvolume files does not match the number of wgmm files');
        for fi = 1:numel(atlFiles)
            subvolReconFlexible(wgmmFiles{fi}, atlFiles{fi}, iniReconFile, atlSubjVol, atlSubjMask, interpMatData, outFiles{fi});
        end
    end
    
    % load subvolume location
    q = load(atlSubvols, 'nfo');
    atlReconLoc = q.nfo.subvolLoc;
    atlSubvolSize = q.nfo.subvolSize;
    
    % get other ingo.
    ini = ini2struct(iniReconFile);
    atlPatchSize = ini.atlPatchSize;
    reconModel = ini.reconModel;
    
    % load volumes
    q = load(wgmmFile, 'wg');
    wg = q.wg;
    wg.params.sigma = wg.wv2sigma;
    
    % get rotation data for this subvolume. do NOT pass in an output file -- we just want this
    % returned. This loads a signle volume, despite the plural (e.g. subvols) in the variable names.
    prepG = strcmp(ini.estep, 'inatlas');
    [a, b, c, d, rangeMins] = vol2subvolInterpData(interpMatData, atlSubjVol, ...
        atlSubjMask, atlReconLoc, atlSubvolSize, '', prepG);
    interpData.subVols = {a};
    interpData.subVolMasks = {b}; 
    interpData.atl2SubjInterpMat = {c};
    interpData.subj2AtlInterpMat = {d};
    
    % get R matrices.
    [interpDataPatches, availableData] = interpData2patchRotData(interpData, atlSubvolSize, 0, atlPatchSize, 1, prepG );
    
    % clean up non-completed patches
    if not(all(availableData(:)))
        fprintf('Using %d available patches\n', sum(availableData(:)));
        assert(any(availableData(:)), 'Could not find available patch');
        interpDataPatches.lenR = interpDataPatches.lenR(availableData, :);
        if prepG, interpDataPatches.lenG = interpDataPatches.lenG(availableData, :); end
        interpDataPatches.yorigAll = interpDataPatches.yorigAll(availableData, :);
        interpDataPatches.yrotmasksAll = interpDataPatches.yrotmasksAll(availableData, :);
        interpDataPatches.ydsmasksAll = interpDataPatches.ydsmasksAll(availableData, :);
        interpDataPatches.boundingBoxes = interpDataPatches.boundingBoxes(availableData, :);
        interpDataPatches.completed = interpDataPatches.completed(availableData, :);
        interpDataPatches.rWeightAll = interpDataPatches.rWeightAll(availableData, :);
        interpDataPatches.ydsmasksAllFullVoxels = interpDataPatches.ydsmasksAllFullVoxels(availableData, :);
    end

    % reconstruct patches
    nanmaskvox2volfn = @(a,b) maskvox2vol(a,b, @nan);
    nPatches = numel(interpDataPatches.yorigAll);
    volRotData = paffine.prepareWgmmRData(interpDataPatches, 1:nPatches, prepG);
    volRotData.K = numel(wg.params.pi);
    if prepG
        Y = zeros(nPatches, size(volRotData.G.data, 1));
        W = zeros(nPatches, size(volRotData.G.data, 1));
        for i = 1:nPatches
            Y(i,:) = volRotData.Y{i} * volRotData.G.data(:, volRotData.G.idx{i})';
            W(i,:) = volRotData.ydsmasks{i} * volRotData.G.data(:, volRotData.G.idx{i})';
        end
        volRotData.extradata = struct('Y', Y, 'W', double(W > ini.estepwthr), 'wts', W, 'K', volRotData.K);
    end
    reconPatches = wg.recon(volRotData, reconModel); %'latentMissingR');
%     reconPatches = cellfunc(@(x) x(:), reconPatches);
    
    % reconstruct subvolume
    rotmasks = interpDataPatches.yrotmasksAll(:);
    patches = cellfunc(nanmaskvox2volfn, reconPatches(:), rotmasks);
    dsmasks = cellfunc(@maskvox2vol, interpDataPatches.ydsmasksAll(:), rotmasks);
    recmasks = cellfunc(@maskvox2vol, volRotData.rWeight, rotmasks);
    recmasks = cellfunc(@(rm, dm) nanmax(rm, dm), recmasks, dsmasks);
    for pi = 1:numel(patches), recmasks{pi}(isnan(patches{pi})) = 0; end
    volStarts = cellfunc(@(x) x(1:3), interpDataPatches.boundingBoxes);  
    
    % determine what region of the subvolume is quilted using the reconPatches
        
%     reconLocsEnd = cellfunc(@(x, y) x  + size(y), ?, patches); 
%     minSubvolLoc = min(cat(1, ?{:}), [], 1);
    minSubvolLoc = rangeMins;
%     maxSubvolLoc = max(cat(1, reconLocsEnd{:}), [], 1);
    patchSizes = cellfunc(@(y) [size(y, 1), size(y, 2), size(y, 3)], patches); % force 3d
    relReconLocsEnd = cellfunc(@(x, y) x  + y, volStarts, patchSizes); 
    subvolAtlSize = max(cat(1, relReconLocsEnd{:}));
    
    % determine the size of the subvolume that will be quilted
%     subvolAtlSize = maxSubvolLoc - minSubvolLoc + 1; 
    
    [reconVol, reconWeight] = patchlib.quiltIrregularPatches(volStarts, patches, ...
        'volSize', subvolAtlSize, 'weightPatches', recmasks);

    reconLoc = minSubvolLoc; 
    
    % save reconPatches.
    save(subjoutFile, 'reconVol', 'reconLoc', 'reconWeight');
