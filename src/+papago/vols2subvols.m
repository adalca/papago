function varargout = vols2subvols(method, data, volnames, atlVolSize, subvolSize, gridSpacing, varargin)
% vols2subvols break down volumes to a set of (potentially overlapping) subvolume columns
%
% subvols = vols2subvols(method, data, atlVolSize, subvolSize, gridSpacing) breaks down volume
% datasets set into subvolume columns. atlVolSize is the size of the common atlas frame. subvolSize
% is the size of the subvolume patch. gridSpacing is the spacing on which to put the subvolumes.
% Potential ways to load volumes includes all the ways used in subspacetools.loadpatches. Notably:
%   method = 'volumes', data = cell{N, K} of volumes.
%   method = 'matfiles', data = cell{N, 1} of matfile pointers
%   volnames = cell{1, K} of volume names. if method is matfile, then this has to be the names in
%       the matfile. Otherwise, This will only act as the names saved to the savefile.
% Note that loading all the subvolume columns in memory can be costly, especially if they are
% overlapping or you have a lot of volumes. We recommend saving the subvolumes and not asking for an
% output argument (see next option).
%
% vols2subvols(method, data, atlVolSize, subvolSize, gridSpacing, savefile) allows for the
% specification of a savefile for each subvolume. the savefile must include absolute path and a '%d'
% to allow for indexing, e.g.: '/path/fname_%d.mat'. the matfile will have a 'subvolumes' cell of
% {NxK} subvolumes saved.
%
% vols2subvols(..., selidx) decides which gridpoints to run. This is useful in several instances.
%   - running in parallel with matfiles. but careful: this might kill the disks!
%   - re-running (or continue running) specific indices.
%
% subvols = vols2subvols(method, data, atlVolSize, subvolSize, gridSpacing)   

    % input parsing
    narginchk(6, 8);

    % setup main gmm grid
    assert(all(gridSpacing > 0));
    subvolOverlap = subvolSize - gridSpacing;
    [subgridlocs, subGridLength] = papago.prepgrid(atlVolSize, subvolSize, subvolOverlap);

    % parse inputs. Need to do this here to have a decent default for idx.
    p = inputParser();
    p.addOptional('savefile', '', @ischar);
    p.addOptional('idx', 1:subGridLength, @isnumeric);
    p.parse(varargin{:});
    savefile = p.Results.savefile;
    idx = p.Results.idx;
    nDims = numel(atlVolSize);
    
    % compute fullgrid.
    fullgridlocs = papago.prepgrid(atlVolSize, ones(1, nDims), 'discrete');

    % get the mapping of the full grid to the closest subgrid point
    centroidlocs = bsxfun(@plus, subgridlocs, floor(gridSpacing/2));
    full2subloc = full2subgridmap(centroidlocs, fullgridlocs, atlVolSize);
    
    % visualize
    % r = randperm(numel(unique(full2subloc)));
    % view3Dopt(reshape(r(full2subloc), atlVolSize)); drawnow;
    
    % prepare output
    if nargout > 0
        varargout{1} = cell(1, subGridLength);
    end
    
    % go through each index, and extract subvolume
    vi = verboseIter(idx, 2);
    while vi.hasNext()
        i = vi.next();
        gridloc = subgridlocs(i, :);

        % extract subvolumes
        % here we're treating the subvolumes as a patch so we can use loadpatches(...)
        % we need to crop the volsize, otherwise we won't get any "patch"
        cropAmount = max((gridloc + subvolSize - 1) - atlVolSize, 0);
        localsubvolSize = subvolSize - cropAmount;
        
        % problem: loadpatches eventually calls subvols2patchvol, which takes the already computed
        % subvolumes and calls patchlib, which takes a long time to do nothing (in this case),
        % essentially. We'll write special cases for 'volumes' and 'matfiles'. TODO: make some other
        % function for this part?
        switch method
            case 'volumes'
                vols = data;
                assert(all(all(cellfun(@(v) all(size(v) == atlVolSize), vols))));
                croprange = arrayfunc(@(loc, sz) loc:(loc+sz-1), gridloc, localsubvolSize);
                subvolumes = cellfunc(@(x) x(croprange{:}), vols);
                
            case 'matfiles'
                matfiles = data;
                croprange = cellfunc(@(loc, sz) loc:sz, gridloc, localsubvolSize);
                subvolumes = subspacetools.matfiles2volumes(matfiles, volnames, croprange);
                
            otherwise               
                subvolumes = subspacetools.loadpatches(method, data{:}, ...
                    localsubvolSize, gridloc, zeros(1, nDims));
                subvolumes = ifelse(iscell(subvolumes), subvolumes, {subvolumes});

                % subvolumes is not a 1xK cell, and each entry is a nVols x prod(localsubvolSize)
                % change it to be a {nVols x K} cell of localsubvolSize sizes
                assert(size(subvolumes{1}, 2) == prod(localsubvolSize));
                subvolumes = cellfunc(@(x) reshape(x', [localsubvolSize, size(x, 1)]), subvolumes);
                subvolumes = cellfunc(@(x) dimsplit(nDims+1, x), subvolumes);
                subvolumes = cellfunc(@(x) x(:), subvolumes);
                subvolumes = [subvolumes{:}];
        end
        assert(all(size(subvolumes{1} > 0)));
        
        % record location information
        nfo.gridloc = gridloc;
        nfo.fulloc = find(full2subloc == i);
        nfo.params = structrich(atlVolSize, subvolSize, gridSpacing, savefile);
        
        % save volumes and information
        if exist('savefile', 'var') && ~isempty('savefile');
            savedata = struct();
            for vi = 1:numel(volnames)
                savedata.(volnames{vi}) = subvolumes(:, vi);
            end
            savedata.nfo = nfo; %#ok<STRNU>
            
            fname = sprintf(savefile, i);
            save(fname, '-struct', 'savedata');
        end

        % write outputs
        if nargout > 0 % watch out, this can be large...
            if i == 1
                q = whos('subvolumes');
                b = q.bytes;
                hrb = humanReadableBytes(b);
                sys.warnif(b > 10*1024^3, 'returning these subvolumes can take %s\n', hrb);
            end
            varargout{1}{i} = subvolumes;
        end
        
        if nargout > 1
            varargout{2}{i} = nfo;
        end
    end
    vi.close();
end

function full2subloc = full2subgridmap(subgridlocs, fullgridlocs, atlVolSize)
% optimal implementation can probably use knowledge of the grids and locations, rather than
% searching for the closest subgrid point for each full grid point, essentially. use gridSpacing,
% etc from parent function.
%
% TODO: move to patchlib, and pass in just gridSpacing and origVolSize to do it.

    % Method 1. This is very slow.
    % full2subloc = zeros(fullGridLength, 1);
    % for i = 1:fullGridLength
    %     [~, full2subloc(i)] = min(pdist2(fullgridlocs(i, :), subgridlocs));
    % end
    
    % Method 2. 
    % Faster, but probably still not optimal.
    % Also, needs some cleaning.
    subgrididx = subvec2ind(atlVolSize, subgridlocs);
    dotsvol = zeros(atlVolSize); 
    dotsvol(subgrididx) = 1;
    [~, closestdotvol] = bwdist(dotsvol);
    
    dotlabels = zeros(atlVolSize); 
    dotlabels(subgrididx) = 1:numel(subgrididx);
    idx = subvec2ind(atlVolSize, fullgridlocs);
    q = closestdotvol(idx); 
    full2subloc = dotlabels(q);
    
    assert(size(full2subloc, 1) == size(fullgridlocs, 1));
end
