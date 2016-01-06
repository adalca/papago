function md2subvols(trainmdfile, mods, volNames, patchSize, gridSpacing, atlVolSize, savefile, mdmatfile)
% MD2SUBVOLS split medical dataset into subvolume columns.
%
% md2subvols(trainmdfile, mods, volNames, patchSize, gridSpacing, atlVolSize, savefile) split
% medical dataset into patch columns. There are two main methods: loading in all volumes (so you
% need to be able to fit all volumes in memory) and chopping the volumes, or loading the right crop
% via matfiles. The latter is achieved by passing in the matfile-modality in the md structure (see
% below). Inputs:
%   trainmdfile e.g. /path/to/md/file
%   mods e.g. "brainIso2Ds%dUs%dsizeReg, brainDs5Us5Reg, brainDs5Us5RegMask" (or maybe all?)
%   volNames e.g. same as mods, or "iso, ds, mask"
%   patchSize e.g. [9, 9, 9];
%   gridSpacing e.g. [5, 5, 5];
%   atlVolSize e.g. [100, 100, 100];
%   savefile e.g. savefile = strrep(fullfile(SYNTHESIS_DATA_PATH, 'adni', 'subvols', 'subvol_%d.mat'), '\', '/');
%
% md2subvols(trainmdfile, mods, volNames, patchSize, gridSpacing, atlVolSize, savefile, mdmatfile)
% allows the specification of a matfile modality
%
% (MCC-ready) numerical inputs can be passed in as strings.
%
% TODO: allow for a mode that just loads in the volume, crops, and stores, one by one (i.e. not
% matfile), and allow to SGE this. This will hit the disks a bit too strong, but might be doable?
%
% See Also papago.vols2subvols, subspacetools.loadpatches

    % parse inputs
    narginchk(7, 8);
    patchSize = makenum(patchSize);
    gridSpacing = makenum(gridSpacing);
    atlVolSize = makenum(atlVolSize);
    if ~iscell(mods), mods = cellfunc(@strtrim, str2cell(mods, ',')); end
    if ~iscell(volNames), volNames = cellfunc(@strtrim, str2cell(volNames, ',')); end

    % load all volumes or matfilepointers
    if nargin == 8
        matfiles = subspacetools.loadpatches('md', trainmdfile, mdmatfile, mods);
        method = {'matfiles', matfiles};
        
    else
        vols = subspacetools.loadpatches('md', trainmdfile, mods);
        method = {'volumes', vols};
    end
    
    % prepare grid
    subvolSize = patchSize + gridSpacing;
    papago.vols2subvols(method{:}, volNames, atlVolSize, subvolSize, gridSpacing, savefile);
end

function x = makenum(x)
    x = ifelse(ischar(x), 'str2num(x)', 'x', true); %#ok<ST2NM>
end
