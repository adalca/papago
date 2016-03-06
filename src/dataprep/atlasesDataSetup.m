%% setup atlases
padAmts = [-1, 0, 10, 30];
names = {'wholevol', 'brain', 'brain_pad10', 'brain_pad30'};

for pi = 1:numel(padAmts)
    p = padAmts(pi);
    name = names{pi};
    
    [BUCKNER_ATLAS_ORIG, bproc] = preprocAtlasPaths('buckner', name, GENERAL_DATA_PATH, ...
        SYNTHESIS_DATA_PATH, 2:7);
    eval(sprintf('BUCKNER_ATLAS_PROC_%s = bproc;', upper(name)));

    [STROKE_ATLAS_ORIG, bproc] = preprocAtlasPaths('stroke', name, GENERAL_DATA_PATH, ...
        SYNTHESIS_DATA_PATH, 7);
    eval(sprintf('STROKE_ATLAS_PROC_%s = bproc;', upper(name)));
end
