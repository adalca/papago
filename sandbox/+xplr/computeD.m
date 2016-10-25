% initialize
setup

%% parameters

atlPatchSize = ones(1, 3) * 9;
atlLoc = LOC_VENTRICLE_EDGE; %LOC_VENTRICLE_EDGE; % LOC_LEFT_CORTEX+10; %LOC_VENTRICLE_EDGE; %LOC_LEFT_CORTEX;
patchColPad = ones(1, 3) * 1;

% train and test datasets
traindataset = 'adni';

% modalities
ds = 5;
us = 5;
isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', ds, us);
dsSubjMod = sprintf('brainDs%dUs%d', ds, us);
dsSubjMaskMod = sprintf('brainDs%dUs%dMask', ds, us);


%% load buckner volumes and prepare volume data
% load ADNI full-subject, and buckner full-dataset column.

% load various data columns from training set
fnames = fullfile(SYNTHESIS_DATA_PATH, traindataset, 'md', [sys.usrname, '_restor_md_*']);
trainmd = loadmd(fnames);
[bucknerIsoPatchCol, layeridx, volidx] = ...
    subspacetools.md2patchcol(trainmd, isoSubjMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjMod, atlPatchSize, atlLoc, patchColPad);
[bucknerDsMaskPatchCol, ~, ~] = ...
    subspacetools.md2patchcol(trainmd, dsSubjMaskMod, atlPatchSize, atlLoc, patchColPad);


%%

% get unique masks
[uniquemask, ia, ic] = unique(bucknerDsMaskPatchCol, 'rows'); 

for i=1:max(ic)
    
    % average the correlation for a given mask
    valid = logical(ic==i);
    diff = bucknerIsoPatchCol(valid,:) - bucknerDsPatchCol(valid,:); 
    corr(:,:,i) = diff'*diff ./ sum(valid); 
    
end
    