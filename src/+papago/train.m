function train(infile, volnames, outfile, patchSize, K, nSamples, model, varargin)
% TRAIN train (weighted) gaussian mixture models 
%
% train(infile, volnames, outfile, patchSize, K, nSamples, model) trains wgmm on the cluster.
% This is meant as as wrapper around papago.learngmm. It sets some default parameters and settings,
% which can mostly be over-written, and requires the most important parameters to be passed in as
% follows
%
%   infile - inputs volumes file, which should have (sub) volumes inside. input files should ahve
%   volumes as well as a 'nfo' field with the data pertaining to how those volumes were created.
%
%   volnames - the volume names inside the matfile necessary here. Usually it's just the data (e.g.
%   isovol) or data and mask (e.g. "dsvol,maskvol")
%   
%   oufile - inputs volumes file, which should have (sub) volumes inside
%   
%   patchSize - 1-by-nDims vector
%   
%   K - scalar, number of clusters 
%   
%   nSampled - number of patches fo be sampled
%       0 - all, (0, 1] - fraction of patches, int > 1 - specific number.
%   
%   model - 'model0' (non-weighted) or 'model3' (latest weighted model)
%
% Allows passing in of numerical parameters as strings, and will convert them to numerical. This
% should help with MCC'ed version of this file where necessary.
%
% We also allow for param/value pairs to overwrite any of the parameters, papago options or wgmm
% options. See papago.learngmm() for the available options. For example, one could pass
% 'Replicates', 5 to force 5 replicates when fitting wgmm. See the parseInputs() function below to
% see all the defaults.
%
% Warning: location has to be relative to this stack of subvolumes
%
%#ok<*ST2NM>

% Main idea: We have some default parameters

    % parameter handling
    [params, opt, wopt] = ...
        parseInputs(infile, volnames, outfile, patchSize, K, nSamples, model, varargin{:});
    
    % load files
    subvols = cellfunc(@(v) load(infile, v), volnames);
    subvols = [subvols{:}];

    % learn gmm.
    tic;
    gmm = papago.learngmm(subvols, params, opt, wopt);
    gmmtoc = toc; %#ok<*NASGU>
    
    % save (bare minimum) gmm
    gmm.X = [];
    gmm.W = [];
    gmm.sigmainv = [];
    save(outfile, 'gmm', 'gmmtoc', 'params', 'opt', 'wgmmopt');
    displ('completed gmm');
end

function [params, opt, wopt] = ...
    parseInputs(infile, volnames, outfile, patchSize, K, nSamples, model, varargin)

    % simple checks and transformations. 
    assert(iseven(numel(varargin)));
    % Note: this requires all Values in Param/Values to have to be numeric.
    varargin(2:2:end) = cellfunc(@(x) makenum(x), varargin(2:2:end));

    % check input files
    assert(sys.isfile(infile));
    volnames = cellfunc(@strtrim, str2cell([volnames, ','], ','));
    w = cellfunc(whos(infile, x), volnames);
    volsfound = sum(cellfun(@(x) numel(x), w));
    assert(volsfound == numel(volnames), 'Did not find all given volumes in input file');
    nVols = prod(w{1}.size);
    assert(ischar(outfile));
    
    % parameters
    patchSize = makenum(patchSize);  
    K = makenum(K);
    nSamples = makenum(nSamples); 
    p = inputParser();
    p.addRequired('patchSize', @(x) isvector(x));
    p.addRequired('K', @isscalar);
    p.addRequired('nSamples', @isscalar);
    p.addParameter('minWeight', 1e-4, @isscalar);
    p.KeepUnmatched = true;
    p.parse(patchSize, K, nSamples, varargin{:});
    params = p.Results;
    params.minWeight = makenum(params.minWeight);

    % papago options
    p = inputParser();
    p.addParameter('subtractpatchmean', true, @islogical);
    p.addParameter('globalinit', false, @islogical);
    p.addParameter('debugisogmm', false, @islogical);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    opt = p.Results;
        
    % wgmm options
    p = inputParser();
    p.addRequired('updateMethod', @(x) ischar(x) && ismember(x, {'model0', 'model3'}));
    p.addParameter('replicates', 3, @isscalar);
    p.addParameter('TolFun', 0.001, @isscalar);
    p.addParameter('regularizationValue', 1e-7, @isscalar);
    p.addParameter('regularizationWeight', nSamples * 5/nVols, @isscalar);
    p.KeepUnmatched = true;
    p.parse(model, varargin{:});
    wopt = p.Results;

end

function x = makenum(x)
    x = ifelse(ischar(x), str2num(x), x);
end
