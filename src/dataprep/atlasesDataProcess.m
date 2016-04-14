%% process BUCKNER atlas 
dsAmt = 2:9;
intNorm = 255;
padAmts = [-1, 0, 10, 30];
names = {'wholevol', 'brain', 'brain_pad10', 'brain_pad30'};

for pi = 1:numel(padAmts)
    p = padAmts(pi);
    fprintf('Running Buckner pad %d\n', p);
    name = names{pi};
    
    [BUCKNER_ATLAS_ORIG, bproc] = preprocAtlas('buckner', GENERAL_DATA_PATH, ...
        SYNTHESIS_DATA_PATH, dsAmt, intNorm, name, p);
    eval(sprintf('BUCKNER_ATLAS_PROC_%s = bproc;', upper(name)));
end


%% process STROKE atlas without pad
dsAmt = 7;
intNorm = 550;

for pi = 1:numel(padAmts)
    p = padAmts(pi);
    fprintf('Running Stroke pad %d\n', p);
    name = names{pi};
    
    [STROKE_ATLAS_ORIG, bproc] = preprocAtlas('stroke', GENERAL_DATA_PATH, ...
        SYNTHESIS_DATA_PATH, dsAmt, intNorm, name, p);
    eval(sprintf('STROKE_ATLAS_PROC_%s = bproc;', upper(name)));
end
