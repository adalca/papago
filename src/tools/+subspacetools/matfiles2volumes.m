function vols = matfiles2volumes(mfptrs, volnames, croprange)
% MATFILES2VOLUMES load several volumes in several matfiles and return in a cell
%   vols = matfiles2volumes(mfptrs, volnames) mfptrs is a size Nx1 cell array of matfile pointers
%   volnames is a size Mx1 cell array of volume names found in each mat file. vols is a size NxM
%   cell of volumes
%
%   vols = matfiles2volumes(mfptrs, volnames, croprange) allows the cropping of each volume as it
%   loads. croprange is a cell such that for a volume, cropvol = volume(croprange{:});

    narginchk(2, 3);
    mfptrs = ifelse(iscell(mfptrs), mfptrs, {mfptrs});
    volnames = ifelse(iscell(volnames), volnames, {volnames});
    
    % cell of matfile pointers to volumes
    vols = cell(numel(mfptrs), numel(volnames));
    for i = 1:numel(mfptrs)
        for v = 1:numel(volnames);
            vol = mfptrs{i}.(volnames{v});
            if nargin == 2
                vols{i, v} = vol;
            else
                vols{i, v} = vol(croprange{:});
            end
        end
    end
    