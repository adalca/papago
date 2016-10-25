function mccReconVolume(fpath, reconnii, dsnii)

    % load template nii
    nii = loadNii(dsnii);
    
    % gather data
    tic
    z3 = wholevol.quiltFromSubvols([fpath, '*.mat']);
    fprintf('done quilting %f\n', toc);
    
    % cut back if necessary
    msize = min(size(nii.img), size(z3));
    z3 = z3(1:msize(1), 1:msize(2), 1:msize(3));
    
    % build up if necessary 
    maxsize = size(nii.img);
    z3(maxsize(1), maxsize(2), maxsize(3)) = 0;
    
    % save
    nii.img = z3;
    saveNii(nii, reconnii);
    