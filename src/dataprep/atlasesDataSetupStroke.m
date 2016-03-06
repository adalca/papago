%% atlasDataSetup
% basic setup of atlas paths
% Requires: GENERAL_DATA_PATH, SYNTHESIS_DATA_PATH
dsAmounts = 2:7;

% stroke paths - general
STROKE_PATH = fullfile(GENERAL_DATA_PATH, 'stroke');
STROKE_ATLAS = fullfile(STROKE_PATH, 'atlases', 'stroke61.nii.gz');
STROKE_ATLAS_SEG = fullfile(STROKE_PATH, 'atlases', 'stroke61_seg.nii.gz');
STROKE_ATLAS_BRAIN = fullfile(STROKE_PATH, 'atlases', 'stroke61_brain.nii.gz');

% stroke paths - processed
STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC = fullfile(SYNTHESIS_DATA_PATH, 'stroke/atlases', 'stroke61_brain_proc.nii.gz');
STROKE_ATLAS_MODS.STROKE_ATLAS_SEG_PROC = fullfile(SYNTHESIS_DATA_PATH, 'stroke/atlases', 'stroke61_seg_proc.nii.gz');

%% DS-US paths
for s = dsAmounts % downsample amount        
    for u = 1:s % upsample amount
        varname = sprintf('STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC_DS%d_US%d', s, u);
        fname = sprintf('stroke61_brain_proc_ds%d_us%d.nii.gz', s, u);
        fullfname = fullfile(SYNTHESIS_DATA_PATH, 'stroke/atlases', fname);
        eval([varname, ' = fullfname;']);
    end
end
