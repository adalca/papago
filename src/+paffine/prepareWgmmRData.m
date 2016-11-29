function q = prepareWgmmRData(interpDataPatches, dataridx)


    % get indices from lengths 
    ridx = lenMatrix2idx(interpDataPatches.lenR);
    assert(~any(cellfun(@isempty, ridx(dataridx))));
    allridx = cat(2, ridx{dataridx});
    
    q.R.data = interpDataPatches.allR(allridx, :);
    gidx = lenMatrix2idx(interpDataPatches.lenG);
    q.G.data = interpDataPatches.allG(:, cat(2, gidx{dataridx}));
    q.R.idx = lenMatrix2idx(interpDataPatches.lenR(dataridx));
    q.G.idx = lenMatrix2idx(interpDataPatches.lenG(dataridx));
    
    q.Y = interpDataPatches.yorigAll(dataridx);
    q.ydsmasks = interpDataPatches.ydsmasksAll(dataridx);
    q.ydsmasksFullVoxels = interpDataPatches.ydsmasksAllFullVoxels(dataridx);
    q.yrotmasks = interpDataPatches.yrotmasksAll(dataridx);

    if isfield(interpDataPatches, 'rWeightAll')
        q.rWeight = interpDataPatches.rWeightAll(dataridx);
    end