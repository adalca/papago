function mccSubvolRecon(gmmFile, subvolFile, iniReconFile, ...
    dsSubjInAtlFile, dsSubjInAtlMaskFile, dsSubjFile, dsSubjWeightFile, ...
    subjCorrFile, subjoutFile)
% for one subject, reconstruct all patches in a subvolume using a gmm
% and combine the patches into a "reconstructed" subvolume via averages,
% along with the weights of the number of patches that voted at each location.
%
% currently, doing simple reconstruction via averages (no mrf)
% currently, require correlation subject file, although tform could also be an opt

    % load 'gmm'
    q = load(gmmFile);
    if isfield(q, 'gmm'), gmm = q.gmm; assert(~isfield(q, 'wg'))
    elseif isfield(q, 'wg'), gmm = q.wg;
    end

    % load subvolume 'nfo'
    load(subvolFile, 'nfo');
    subvolLoc = nfo.subvolLoc;
    subvolSize = nfo.subvolSize;

    % load method from some sort of settings (?)
    ini = ini2struct(iniReconFile);
    crmethod = ini.dirn; % 'inverse';
    nPatchReconPerLoc = ini.mrf.nPatchReconPerLoc; % 1;
    atlPatchSize = ini.atlPatchSize;

    for i = 1:size(gmm.mu, 1)
        gmm.sigma(:,:,i) = gmm.sigma(:,:,i) + eye(size(gmm.mu,2)) * ini.regval;
    end
    
    % load volumes
    dsSubjInAtlNii = loadNii(dsSubjInAtlFile);
    dsSubjInAtlVol = nii2vol(dsSubjInAtlNii); dsSubjInAtlNii.img = [];
    dsSubjInAtlMaskVol = nii2vol(loadNii(dsSubjInAtlMaskFile));
    dsSubjNii = loadNii(dsSubjFile);
    dsSubjVol = nii2vol(dsSubjNii); dsSubjNii.img = [];
    dsSubjWeightVol = nii2vol(loadNii(dsSubjWeightFile));

    % check which type of file is passed: a corr or a tform
    istform = numel(whos(matfile(subjCorrFile), 'tform')) > 0;
    if istform
        subjInAtlTform  = load(subjCorrFile, 'tform');
       
        subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
        atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
        tform = subjInAtlTform.tform;
        atlVolSize = size(dsSubjInAtlVol);
        subjLoc2AtlSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims);
        atlLoc2SubjSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims, 'backward');

        % reconstruct
        [reconVol, reconLoc, reconWeight] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, ...
            crmethod, nPatchReconPerLoc, dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubjVol, ...
            dsSubjWeightVol, atlLoc2SubjSpace, subjLoc2AtlSpace); %#ok<ASGLU>
        
    else
        load(subjCorrFile, 'atlLoc2SubjSpace', 'subjLoc2AtlSpace', 'atl2subjR');

        % reconstruct
        [reconVol, reconLoc, reconWeight] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, ...
            crmethod, nPatchReconPerLoc, dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubjVol, ...
            dsSubjWeightVol, atlLoc2SubjSpace, subjLoc2AtlSpace, atl2subjR); %#ok<ASGLU>
    end

    % save reconPatches.
    save(subjoutFile, 'reconVol', 'reconLoc', 'reconWeight');
