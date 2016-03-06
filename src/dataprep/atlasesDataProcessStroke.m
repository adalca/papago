dsAmounts = 2:7;
intensityNorm = 550;

%% process atlas
normalizeNii(STROKE_ATLAS_BRAIN, STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC, intensityNorm);
padNii(STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC, STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC, [tempname, '.nii.gz'], 10);
[croppedVol, cropMask, cropArray, bBox] = boundingBoxNii(STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC, STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC);

% crop segmentation file
segnii = loadNii(STROKE_ATLAS_SEG);

segcrop = segnii.img(cropArray{:});
niinew = make_nii(segcrop);
niinew.hdr.dime.pixdim(2:4) = segnii.hdr.dime.pixdim(2:4);
saveNii(niinew, STROKE_ATLAS_MODS.STROKE_ATLAS_SEG_PROC);
padNii(STROKE_ATLAS_MODS.STROKE_ATLAS_SEG_PROC, STROKE_ATLAS_MODS.STROKE_ATLAS_SEG_PROC, [tempname, '.nii.gz'], 10);


%% resample (all) atlases to match downsample + upsample process 
atlnii = loadNii(STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC);
atlvol = atlnii.img;

for s = dsAmounts % downsample amount        
    for u = 1:s % upsample amount
        
        % downsample atlas nii
        sz = round(size(atlvol) ./ s * u);
        cmd = 'vol = volresize(atlvol, sz);';
        disp(cmd);
        eval(cmd);
        sunii = make_nii(vol);
        sunii.hdr.dime.pixdim(2:4) = atlnii.hdr.dime.pixdim(2:4) ./ u * s;
        
        % save
        varname = sprintf('STROKE_ATLAS_MODS.STROKE_ATLAS_BRAIN_PROC_DS%d_US%d', s, u);
        cmd = sprintf('saveNii(sunii, %s)', varname);
        disp(cmd);
        eval(cmd);
    end
end
