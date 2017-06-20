function params = mstep(wg, data)
% m-step (parameter updates)
    
    % parse inputs
	narginchk(2, 2);
    switch wg.opts.model.name
        case 'model0'
            params = wg.mstepModel0(data);
        case 'model1'
            params = wg.mstepModel1(data);
        case 'model3'
            params = wg.mstepModel3(data);
        case 'model4exp'
            params = wg.mstepModel4exp(data);
        case 'model5'
            params = wg.mstepModel5(data);
        case 'latentSubspace'
            params = wg.mstepLatentSubspace(data);
        case 'gpuLatentSubspace'
            params = wg.mstepGPULatentSubspace(data);
        case 'latentSubspaceR'
            tdata = struct();
            tdata.Y = data.atlY;
            tdata.W = data.atlW;
            tdata.K = data.K;
            params = wg.mstepLatentSubspace(tdata);
        case 'wLatentSubspace'
            params = wg.mstepWLatentSubspace(data);
        case 'latentMissing'
            params = wg.mstepLatentMissing(data);
        case 'latentMissingR'
            params = wg.mstepLatentMissingR(data);
    end
    
    % check every element
    for f = fieldnames(params)'
        assert(isclean(params.(f{:})));
    end
end