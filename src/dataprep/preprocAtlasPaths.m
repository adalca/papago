function [orig, proc] = preprocAtlasPaths(atlType, procFolder, genPath, procPath, dsAmounts)
% Example
% [orig, atlMods] = prepAtlasPaths('buckner', 'atlases', GENERAL_DATA_PATH, SYNTHESIS_DATA_PATH, 2:7)
   
    if ~strcmp(genPath(end), filesep), genPath(end+1) = '/'; end
    if ~strcmp(procPath(end), filesep), procPath(end+1) = '/'; end

    % atlas name
    atlname = sprintf('%s61', atlType);

    % general paths 
    orig.PATH = fullfile(genPath, atlType);
    orig.ATLAS = fullfile(orig.PATH, 'atlases', [atlname, '.nii.gz']);
    orig.SEG = fullfile(orig.PATH, 'atlases', [atlname, '_seg.nii.gz']);
    orig.BRAIN = fullfile(orig.PATH, 'atlases', [atlname, '_brain.nii.gz']);

    % processed paths 
    proc.PATH = sprintf('%s%s/atlases/%s', procPath, atlType, procFolder);
    base = sprintf('%s%s/atlases/%s/%s', procPath, atlType, procFolder, atlname);
    proc.ATLAS = sprintf('%s_proc.nii.gz', base);
    proc.BRAIN = sprintf('%s_brain_proc.nii.gz', base);
    proc.SEG = sprintf('%s_seg_proc.nii.gz', base);

    % processed DsUs paths
    for s = dsAmounts % downsample amount        
        for u = 1:s % upsample amount
            fieldname = sprintf('ATLAS_DS%d_US%d', s, u);
            fname = sprintf('_proc_ds%d_us%d.nii.gz', s, u);
            fullfname = sprintf('%s%s', base, fname);
            proc.(fieldname) = fullfname;
            
            fieldname = sprintf('BRAIN_DS%d_US%d', s, u);
            fname = sprintf('_brain_proc_ds%d_us%d.nii.gz', s, u);
            fullfname = sprintf('%s%s', base, fname);
            proc.(fieldname) = fullfname;
            
            fieldname = sprintf('SEG_DS%d_US%d', s, u);
            fname = sprintf('_seg_proc_ds%d_us%d.nii.gz', s, u);
            fullfname = sprintf('%s%s', base, fname);
            proc.(fieldname) = fullfname;
        end
    end
end

