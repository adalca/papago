function [orig, proc] = preprocAtlas(atlType, genPath, procPath, dsAmounts, intensityNorm, procFolder, padAmount)
% preprocAtlas(atlType, GENERAL_DATA_PATH, SYNTHESIS_DATA_PATH, dsAmounts, intensityNorm)
%
% set padAmount to < 0 (e.g. -1) to avoid croppping at all!

    % get paths
    [orig, proc] = preprocAtlasPaths(atlType, procFolder, genPath, procPath, dsAmounts);
    
    %% process atlas
    assert(all(size(nii2vol(orig.ATLAS)) == size(nii2vol(orig.BRAIN))));
    assert(all(size(nii2vol(orig.ATLAS)) == size(nii2vol(orig.SEG))));
    
    % normalize brain
    normalizeNii(orig.ATLAS, proc.ATLAS, intensityNorm);
    normalizeNii(orig.BRAIN, proc.BRAIN, intensityNorm);
    
    if padAmount >= 0 
        % crop brain
        [~, ~, cropArray, ~] = boundingBoxNii(proc.BRAIN, proc.BRAIN);

        % pad brain
        prePadAmount = 0;
        postPadAmount = 0;
        if padAmount > 0
            volSize = size(nii2vol(proc.ATLAS));
            padNii(proc.BRAIN, proc.BRAIN, [], padAmount);
            p = padAmount;
            prePadAmount = max(padAmount - cellfun(@(x) x(1), cropArray) + 1, 0);
            postPadAmount = max(cellfun(@(x) x(end), cropArray) + padAmount - volSize, 0);
            cropArray = cellfunc(@(c, v) [max(c(1)-p, 1):c(1)-1, c, c(end)+1:min(c(end)+p, v)], ...
                cropArray, mat2cellsplit(volSize));
        end

        % crop segmentation and whole atlas file
        cropMod(proc.ATLAS, cropArray, proc.ATLAS)
        cropMod(orig.SEG, cropArray, proc.SEG)

        % pad files
        padNii(proc.ATLAS, proc.ATLAS, [], prePadAmount, 'pre');
        padNii(proc.ATLAS, proc.ATLAS, [], postPadAmount, 'post');
        padNii(proc.SEG, proc.SEG, [], prePadAmount, 'pre');
        padNii(proc.SEG, proc.SEG, [], postPadAmount, 'post');
    
    else
        copyfile(orig.SEG, proc.SEG)
    end

    %% resample (all) atlases to match downsample + upsample process 
    dsusproc(dsAmounts, proc, 'ATLAS');
    dsusproc(dsAmounts, proc, 'BRAIN');
    dsusproc(dsAmounts, proc, 'SEG');
    
    %% save
    atlMods = proc;
    save(sprintf('%s/%s', proc.PATH, 'atlMods'), 'atlMods');

end

function dsusproc(dsAmounts, proc, type) 
    brainnii = loadNii(proc.(type));
    brainvol = brainnii.img;

    for s = dsAmounts % downsample amount        
        for u = 1:s % upsample amount

            % downsample atlas nii
            sz = round(size(brainvol) ./ s * u);
            vol = volresize(brainvol, sz);
            sunii = make_nii(vol);
            sunii.hdr.dime.pixdim(2:4) = brainnii.hdr.dime.pixdim(2:4) ./ u * s;

            % save
            fieldname = sprintf('%s_DS%d_US%d', type, s, u);
            saveNii(sunii, proc.(fieldname));
        end
    end
end

function cropMod(inMod, cropArray, saveMod)
    segnii = loadNii(inMod);
    segcrop = segnii.img(cropArray{:});
    niinew = make_nii(segcrop);
    niinew.hdr.dime.pixdim(2:4) = segnii.hdr.dime.pixdim(2:4);
    saveNii(niinew, saveMod);
end
