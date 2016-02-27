dsAmounts = 2:7;
intensityNorm = 255;

%% process atlas
normalizeNii(BUCKNER_ATLAS_BRAIN, BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC, intensityNorm);
[croppedVol, cropMask, cropArray, bBox] = boundingBoxNii(BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC, BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC);

% crop segmentation file
segnii = loadNii(BUCKNER_ATLAS_SEG);
segcrop = segnii.img(cropArray{:});
niinew = make_nii(segcrop);
niinew.hdr.dime.pixdim(2:4) = segnii.hdr.dime.pixdim(2:4);
saveNii(niinew, BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_SEG_PROC);

%% resample (all) atlases to match downsample + upsample process 
atlnii = loadNii(BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC);
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
        varname = sprintf('BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', s, u);
        cmd = sprintf('saveNii(sunii, %s)', varname);
        disp(cmd);
        eval(cmd);
    end
end
