function processSegFlairAtl(segfile, maskfile)
% flair is ds
% assume t1file is iso

    % estimate dsrate.
    nii = loadNii(segfile);
    nii.img = double(nii.img > 0);
    saveNii(nii, maskfile);