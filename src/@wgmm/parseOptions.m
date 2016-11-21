function opts = optionDefaults()
% option defaults
    
    opts.model.name = 'model0';
    opts.init.method = 'exemplar';
    opts.replicates = 10;
    opts.maxIter = 10;
    opts.regularizationValue = 1e-7;
    opts.TolFun = 0.01;
    opts.verbose = false;
end
