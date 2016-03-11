function nii2subvols(niifile, volName, patchSize, gridSpacing, atlVolSize, savefile)
% NII2SUBVOLS split medical dataset into subvolume columns.
%
% nii2subvols(niifile, volName, patchSize, gridSpacing, atlVolSize, savefile) split medical nifti
% into patch columns. There are two main methods: loading in volume and chopping the volumes
%
% Inputs:
%   trainmdfile e.g. /path/to/md/file
%   niifile e.g. /path/to/buckner01_iso_2_ds5_us5_size_reg.nii.gz
%   volNames e.g. 'subvolume'
%   patchSize e.g. [9, 9, 9];
%   gridSpacing e.g. [5, 5, 5];
%   atlVolSize e.g. [100, 100, 100];
%   savefile e.g. savefile = strrep(fullfile(SYNTHESIS_DATA_PATH, 'adni', 'subvols', 'subvol_%d.mat'), '\', '/');
%
% md2subvols(trainmdfile, mods, volNames, patchSize, gridSpacing, atlVolSize, savefile, mdmatfile)
% allows the specification of a matfile modality
%
% (MCC-ready) numerical inputs can be passed in as strings. % TODO: change this, make any
% mcc-readyness outside of this function (?).
%
% TODO: allow for a mode that just loads in the volume, crops, and stores, one by one (i.e. not
% matfile), and allow to SGE this. This will hit the disks a bit too strong, but might be doable?
%
% See Also papago.vols2subvols, subspacetools.loadpatches

    % parse inputs
    narginchk(6, 6);
    patchSize = makenum(patchSize);
    gridSpacing = makenum(gridSpacing);
    atlVolSize = makenum(atlVolSize);

    % load volume
    vol = nii2vol(niifile);
    method = {'volumes', {vol}};
    
    % prepare grid
    subvolSize = patchSize + gridSpacing - 1; 
    papago.vols2subvols(method{:}, {volName}, atlVolSize, subvolSize, gridSpacing, savefile);
end

function x = makenum(x)
    x = ifelse(ischar(x), 'str2num(x)', 'x', true); 
end
