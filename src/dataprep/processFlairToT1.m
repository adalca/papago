function processFlairToT1(flairfile, t1file, intensityNorm, flairfileout, flairfilemaskout)
% flair is ds
% assume t1file is iso

    % estimate dsrate.
    flairnii = loadNii(flairfile);
    t1nii = loadNii(t1file);
    dsRate = round(flairnii.hdr.dime.pixdim(4) ./ flairnii.hdr.dime.pixdim(2));
    assert(dsRate == (round(flairnii.hdr.dime.pixdim(4) ./ flairnii.hdr.dime.pixdim(3))));
    
    % us rates
    usRates = 1:dsRate;
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    flairnii.img = double(flairnii.img) ./ intensityNorm;
    [flairusvol, flairusmask] = upsampleNii(flairnii, [], [], 'linear', 0, [1, 1, dsRate], true);
    
    %% register
    movingDims = flairnii.hdr.dime.pixdim(2:4);
    rMoving = imref3d(size(flairusvol), movingDims(2), movingDims(1), movingDims(3));
    fixedDims = t1nii.hdr.dime.pixdim(2:4);
    rFixed = imref3d(size(t1nii.img), fixedDims(2), fixedDims(1), fixedDims(3));
    
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(flairusvol, rMoving, t1nii.img, rFixed, 'rigid', optimizer, metric);
    flairusvol_reg = imwarp(flairusvol, tform, 'OutputView', rFixed);
    flairusmask_reg = imwarp(double(flairusmask), tform, 'OutputView', rFixed);
    
    %% save 
    flairoutnii = t1nii;
    flairoutnii.img = flairusvol_reg;
    saveNii(flairoutnii, flairfileout);
    flairmasknii = t1nii;
    flairmasknii.img = flairusmask_reg;
    saveNii(flairmasknii, flairfilemaskout);
    