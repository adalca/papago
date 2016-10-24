function vol = quiltVolFromSubvols(filenames)
%
% example
% >> datapath = '/data/vision/polina/scratch/adalca/patchSynthesis/data'
% >> b01filenames = fullfile(datapth, 'buckner01/subvolRecons/dec2015/model0/K5/*.mat')
% b01recon = quiltVolFromSubvols(b01filenames)

    d = sys.fulldir(filenames);

    reconVols = cell(1, numel(d)); 
    reconLocs = reconVols; 
    reconWeights = reconVols; 

    vi = verboseIter(1:numel(d), 2);
    while vi.hasNext()
        i = vi.next();
        q = load(d(i).name); 
        reconVols{i} = q.reconVol; 
        reconLocs{i} = q.reconLoc; 
        reconWeights{i} = q.reconWeight; 
    end
    vi.close();

    % quilt from irregular patches
    vol = patchlib.quiltIrregularPatches(reconLocs, reconVols, 'weightPatches', reconWeights); 
end
