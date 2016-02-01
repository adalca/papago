%% atlasDataSetup
% basic setup of atlas paths
% Requires: GENERAL_DATA_PATH, SYNTHESIS_DATA_PATH
dsAmounts = 2:7;

% buckner paths - general
BUCKNER_PATH = fullfile(GENERAL_DATA_PATH, 'buckner');
BUCKNER_ATLAS = fullfile(BUCKNER_PATH, 'atlases', 'buckner61.nii.gz');
BUCKNER_ATLAS_SEG = fullfile(BUCKNER_PATH, 'atlases', 'buckner61_seg.nii.gz');
BUCKNER_ATLAS_BRAIN = fullfile(BUCKNER_PATH, 'atlases', 'buckner61_brain.nii.gz');

% buckner paths - processed
BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC = fullfile(SYNTHESIS_DATA_PATH, 'buckner/atlases', 'buckner61_brain_proc.nii.gz');
BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_SEG_PROC = fullfile(SYNTHESIS_DATA_PATH, 'buckner/atlases', 'buckner61_seg_proc.nii.gz');

%% DS-US paths
for s = dsAmounts % downsample amount        
    for u = 1:s % upsample amount
        varname = sprintf('BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', s, u);
        fname = sprintf('buckner61_brain_proc_ds%d_us%d.nii.gz', s, u);
        fullfname = fullfile(SYNTHESIS_DATA_PATH, 'buckner/atlases', fname);
        eval([varname, ' = fullfname;']);
    end
end
