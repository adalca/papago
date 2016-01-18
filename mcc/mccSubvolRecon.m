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
    load(gmmFile);
    
    % load subvolume 'nfo' 
    load(subvolFile, 'nfo'); 
    subvolLoc = nfo.gridloc;
    subvolSize = nfo.params.subvolSize;

    % load method from some sort of settings (?)
    ini = ini2struct(iniReconFile);
    crmethod = ini.dirn; % 'inverse';
    nPatchReconPerLoc = ini.mrf.nPatchReconPerLoc; % 1;
    atlPatchSize = ini.atlPatchSize;

    % load volumes
    dsSubjInAtlVol = nii2vol(loadNii(dsSubjInAtlFile));
    dsSubjInAtlMaskVol = nii2vol(loadNii(dsSubjInAtlMaskFile));
    dsSubjVol = nii2vol(loadNii(dsSubjFile));
    dsSubjWeightVol = nii2vol(loadNii(dsSubjWeightFile));
    
    load(subjCorrFile, 'atlLoc2SubjSpace', 'subjLoc2AtlSpace');

    % reconstruct
    [reconVol, reconLoc, reconWeight] = papago.subvolRecon(gmm, subvolLoc, subvolSize, atlPatchSize, ...
        crmethod, nPatchReconPerLoc, dsSubjInAtlVol, dsSubjInAtlMaskVol, dsSubjVol, ...
        dsSubjWeightVol, atlLoc2SubjSpace, subjLoc2AtlSpace); %#ok<ASGLU>

    % save reconPatches.
    save(subjoutFile, 'reconVol', 'reconLoc', 'reconWeight');
