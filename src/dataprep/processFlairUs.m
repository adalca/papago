function processFlairUs(flairfile, intensityNorm, flairfileout, flairfilemask)
% flair is ds
% assume t1file is iso

    % estimate dsrate.
    flairnii = loadNii(flairfile);
    dsRate = round(flairnii.hdr.dime.pixdim(4) ./ flairnii.hdr.dime.pixdim(2));
    assert(dsRate == (round(flairnii.hdr.dime.pixdim(4) ./ flairnii.hdr.dime.pixdim(3))));
    
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    flairnii.img = double(flairnii.img) ./ intensityNorm;
    [flairusvol, flairusmask] = upsampleNii(flairnii, [], [], 'linear', 0, [1, 1, dsRate], true);
    
    %% save 
    flairoutnii = make_nii(flairusvol);
    saveNii(flairoutnii, flairfileout);
    flairmasknii = make_nii(flairusmask*1);
    saveNii(flairmasknii, flairfilemask);
    
