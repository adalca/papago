function [vol, countvol] = quiltSubvols(subVolPath, quiltVolumeFilename)

    d = sys.fulldir(subVolPath); 
    fprintf('Found %d subvolumes\n', numel(d));
    
    if numel(d) == 0
        vol = [];
        countvol = [];
    else

        reconVols = cell(1, numel(d)); 
        reconLocs = reconVols; 
        reconWeights = reconVols; 

        for i = 1:numel(d)
            q=load(d(i).name); 
            reconVols{i}=q.reconVol; 
            reconLocs{i}=q.reconLoc; 
            reconWeights{i}=q.reconWeight; 
        end


        [vol, countvol] = patchlib.quiltIrregularPatches(reconLocs, reconVols, 'weightPatches', reconWeights);
    end
    
    if exist('quiltVolumeFilename', 'var')
        save(quiltVolumeFilename, 'vol', 'countvol'); 
    end    
