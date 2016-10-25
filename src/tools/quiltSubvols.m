function quiltSubvols(subVolPath, quiltVolumeFilename)

d = sys.fulldir(subVolPath); 

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

save(quiltVolumeFilename, 'vol', 'countvol'); 
