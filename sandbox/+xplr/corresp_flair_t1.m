%% explore a small subset of ADNI data where we have both flair and t1.
% explore:
% 1. making a FLAIR atlas from T1 data
% 2. Transfering sigmas from one modality to another.

% initialize
setup

%% parameters
atlPatchSize = ones(1, 3) * 9; 
atlLoc = LOC_VENTRICLE_EDGE; 
%atlLoc = [106,96,125]
reconSubj = 1; %1, 3
patchColPad = ones(1, 3) * 2;

% train and test datasets
testdataset = 'adni';

% recon params
subvolLoc = atlLoc - patchColPad;
subvolSize = atlPatchSize + (2 * patchColPad + 1);

% isoSubjInAtlMod = sprintf('Iso2Ds%dUs%dsizeReg', ds, us);
% dsSubjInAtlMod = sprintf('Ds%dUs%dReg', ds, us);
% dsSubjInAtlMaskMod = sprintf('Ds%dUs%dRegMask', ds, us);
% dsSubjInAtlMatMod = sprintf('Ds%dUs%dRegMat', ds, ds); % note: meant to be ds, ds !
% dsInterpSubjInAtlMod = sprintf('Ds%dUs%dInterpReg', ds, us);
% 
% dsSubjMod = sprintf('Ds%dUs%d', ds, us);
% dsSubjMaskMod = sprintf('Ds%dUs%dMask', ds, us);
% isoSubjMod = sprintf('Iso2Ds%dUs%dsize', ds, us);

%% get md
md = medicalDataset();
md.addRequiredModality('flairInAtl', '%s_flair_to__iso_2_ds5_us5_size_reg.nii.gz');
md.addRequiredModality('flairMaskInAtl', '%s_flair_to__iso_2_ds5_us5_size_reg_dsmask.nii.gz');
md.addRequiredModality('t1InAtl', '%s_iso_2_ds5_us5_size_reg.nii.gz');
md.build([SYNTHESIS_DATA_PATH, 'FLAIR_atlas_construction/atlmk'], []);
%%
for i = 1:md.getNumSubjects, 
    nii = md.loadModality('flairInAtl', i); 
    nii.img = double(nii.img); 
    nii.hdr.dime.datatype = 64;nii.hdr.dime.bitpix = 64;
    saveNii(nii, md.getModality('flairInAtl', i)); % an't use saveModality since modality is protected
    nii = md.loadModality('flairMaskInAtl', i); 
    nii.img = double(nii.img); 
    nii.hdr.dime.datatype = 64;nii.hdr.dime.bitpix = 64;
    saveNii(nii, md.getModality('flairMaskInAtl', i)); % an't use saveModality since modality is protected
    nii = md.loadModality('t1InAtl', i); 
    nii.img = double(nii.img); 
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    saveNii(nii, md.getModality('t1InAtl', i)); 
end

%% get atlas (T1)
atlnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/brain_pad10/buckner61_brain_proc.nii.gz']);
atlnii.img./55*40.5; % histogram correction
atlsegnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/brain_pad10/buckner61_seg_proc.nii.gz']);

flatlnii = loadNii([SYNTHESIS_DATA_PATH, 'stroke/atlases/brain_pad10/stroke61_brain_proc.nii.gz']);

%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load various data columns from training set
[flairPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(md, 'flairInAtl', atlPatchSize, atlLoc, patchColPad);
[flairMaskPatchCol, ~, volidx] = ...
    subspacetools.md2patchcol(md, 'flairMaskInAtl', atlPatchSize, atlLoc, patchColPad);
[t1PatchCol, ~, ~] = ...
    subspacetools.md2patchcol(md, 't1InAtl', atlPatchSize, atlLoc, patchColPad);


% % load selected ADNI subject volumes
% fnames = fullfile(SYNTHESIS_DATA_PATH, testdataset, 'md', [sys.usrname, '_restor_md_*']);
% testmd = loadmd(fnames);
% dsSubjNii = testmd.loadModality(dsSubjMod, reconSubj);
% dsSubjVol = double(dsSubjNii.img);
% dsSubjWeightVol = logical(testmd.loadVolume(dsSubjMaskMod, reconSubj));
% isoSubjInAtlNii = testmd.loadModality(isoSubjInAtlMod, reconSubj);
% dsSubjInAtlNii = testmd.loadModality(dsSubjInAtlMod, reconSubj);
% dsSubjInAtlMaskVol = testmd.loadVolume(dsSubjInAtlMaskMod, reconSubj);
% subjInAtlTform = load(testmd.getModality(dsSubjInAtlMatMod, reconSubj));
% 
% % prepare necessary inputs for conditional-based reconstruction
% subjDims = dsSubjNii.hdr.dime.pixdim(2:4);
% atlDims = dsSubjInAtlNii.hdr.dime.pixdim(2:4);
% tform = subjInAtlTform.tform;
% atlVolSize = size(dsSubjInAtlNii.img);
% subjLoc2AtlSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims);
% atlLoc2SubjSpace = tform2cor3d(tform, size(dsSubjVol), subjDims, atlVolSize, atlDims, 'backward');

%%
disp('knn quilt'); 
warning('cleanup result by label')

% crop the subvolumes from the iso, ds and mask volumes
atlSubvol = cropVolume(double(atlnii.img), subvolLoc, subvolLoc + subvolSize - 1);

% divide the volumes into patch libraries
libPatches = patchlib.vol2lib(atlSubvol, atlPatchSize); 

% find the matching patches 
K = 10;
[pIdx, dist] = knnsearch(t1PatchCol, libPatches, 'K', K); 
retrievedPatches = zeros(size(pIdx, 1), size(libPatches, 2), K);
retrievedMaskPatches = zeros(size(pIdx, 1), size(libPatches, 2), K);
for i = 1:K
    retrievedPatches(:,:,i) = flairPatchCol(pIdx(:,1), :); 
    retrievedMaskPatches(:,:,i) = flairMaskPatchCol(pIdx(:,i), :); 
end

% quilt the patches
wt = bsxfun(@times, exp(-permute(dist, [1, 3, 2])), retrievedMaskPatches);
wt = retrievedMaskPatches + 0.01;
reconVolw = patchlib.quilt(retrievedPatches, size(atlSubvol) - atlPatchSize + 1, atlPatchSize, 'nnWeights', wt); 
reconVol = patchlib.quilt(retrievedPatches, size(atlSubvol) - atlPatchSize + 1, atlPatchSize); 

% view the results
% linearlly interpolated subvolume, ground truth iso subvolume, reconstruced Volume
flatlSubvol = cropVolume(double(flatlnii.img), subvolLoc, subvolLoc + subvolSize - 1);
segatlSubvol = cropVolume(double(atlsegnii.img), subvolLoc, subvolLoc + subvolSize - 1);
view3Dopt(atlSubvol, reconVolw, reconVol, flatlSubvol, segatlSubvol)

