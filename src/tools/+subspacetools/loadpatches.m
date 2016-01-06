function varargout = loadpatches(loadmethod, varargin)
% LOADPATCHES load patches in various ways for the subspace project.
% 
% varargout = loadpatches(loadmethod, varargin) load patches (or just prepare the relevant loading
% structure) for the subspace project. The exact function signature depends on the loadmethod. The
% loadmethods can currently be:
%   'niifiles' - load from a cell array of nifti files
%   'niistruct' - load from a nifti struct
%   'md' - load from a medicalDataset object pointing to matfiles or different modalities
%   'matfiles' - load from a cell of mat-filenames or matfiles pointers
%   'volumes' - load from a cell of numeric volumes
% The detailed function signatures are available below each subfunction.
%
% TODOs: 
%   Add (online) volume crop support with all methods (currently just in 'matfiles')
%   Add comments based on comments in each function.
%   perhaps make (some of) this part of medicalDataset project or patchlib project?
%
% Contact: adalca@csail

    nout = nargout; 
    varargout = cell(1, nout);

    % prepare the structure first
    switch loadmethod
        case 'niifiles'
            [varargout{:}] = load_niifiles(varargin{:});

        case 'niistruct'
            [varargout{:}] = load_niistruct(varargin{:});

        case 'md'
            [varargout{:}] = load_md(varargin{:});
            
        case 'matfiles'
            [varargout{:}] = load_matfiles(varargin{:});
            
        case 'volumes'
            [varargout{:}] = load_volumes(varargin{:});
            
        otherwise
            error('loadpatches: Unknown load method.');
    end
end

function varargout = load_niifiles(niifiles, varargin)
% vols = load_niifiles(niifiles{N,K})
% fnvols = load_niifiles(niifiles{N,K}, fnvols)
% patches = load_niifiles(niifiles, patchSize, locations)
% [patches, layeridx, volidx] = load_niifiles(niifiles, patchSize, locations)
% [patches, ...] = load_niifiles(niifiles, patchSize, locations, locpad)
% [patches, layeridx, volidx, vols] = load_niifiles(niifiles, patchSize, locations, ...) 
   
    % input checking
    narginchk(1, 4);
    if nargin == 3, varargin{4} = varargin{2} * 0; end % have locpad
    nargoutchk(1, 4);
    isc = iscell(niifiles);
    niifiles = ifelse(isc, niifiles, {niifiles});
    dopatches = nargin > 2;
    dovolload = ~dopatches || nargout == 4;
    
    % load all volumes if desired
    if dovolload 
        if nargin == 2
            fn = varargin{1};
            loadfn = @(x) fn(nii2vol(loadNii(x)));
        else % nargin is 1 or 3, 4
            loadfn = @(x) nii2vol(loadNii(x));
        end 
        vols = cellfunc(loadfn, niifiles);
    end
        
    % just volumes
    if ~dopatches && dovolload 
        varargout = {vols};
        
    % just patches
    elseif dopatches && ~dovolload
        colcmd = @(nii) subspacetools.nii2patchcol(nii, varargin{:});
        
        patches = cell(1, size(niifiles, 2));
        for k = 1:size(niifiles, 2)
            [patchesc, layeridx, volidx] = cellfunc(colcmd, niifiles);
            patches{k} = cat(1, patchesc{:});
        end 
        patches = ifelse(isc, patches, patches{1}); % if niifiles was a cell, keep patches as cell.
        varargout = {patches, layeridx, volidx};
        varargout = varargout(1:nargout);
       
    % all volumes & patches
    else
        assert(dopatches && dovolload, 'Cannot have no patches and no volume');
        [patches, layeridx, volidx] = load_volumes(nargout, vols, varargin{2:end});
        varargout = {patches, layeridx, volidx, vols};
    end
end

function varargout = load_niistruct(varargin) 
% vols = load_niistruct(niistruct, volnames)  % 3 inputs
% vols = load_niistruct(niifile, structname, volnames)  % 4 inputs
% [patches, layeridx, volidx] = load_niistruct(..., patchSize, locations)  % 5 inputs or 6 inputs
% [patches, layeridx, volidx] = load_niistruct(..., patchSize, locations, locpad) % 6 inputs or 7 inputs
% [patches, layeridx, volidx, vols] = load_niistruct(..., patchSize, locations, ...) 
% 
% forces loading of all volumes no matter what

    assert(ismember(nargin, [2, 3, 5, 6]));
    nargoutchk(1, 4);
    
    if ischar(varargin{1}) % niifile
        [niifile, structname, volnames] = varargin{1:3};
        l = load(niifile);
        niistruct = l.(structname);
        othervars = varargin(4:end);
        
    else % niistruct
        [niistruct, volnames] = varargin{1:2};
        othervars = varargin(3:end);
    end
    
    if numel(othervars) == 2, othervars{3} = othervars{1} * 0; end
    
    % forces loading of all volumes no matter what
    vols = cell(1, numel(volnames));
    for v = 1:numel(volnames)
        [vols{:, v}] = arrayfunc(@(i) niistruct.(volnames){i}.img, 1:length(niistruct.(volnames)))'; 
    end
    
    % prepare output
    varargout = cell(nargout);
    [varargout{:}] = load_volumes(vols, othervars{:});
end

function varargout = load_md(varargin)
% mfptrs = load_md(md, matfile_modality);
% [mfptrs, vols] = load_md(md, matfile_modality, volumes2load); if nargout == 2
% [patches, layeridx, volidx, mfptrs, vols] = load_md(..., patchSize, locations);
% [patches, layeridx, volidx, mfptrs, vols] = load_md(..., patchSize, locations, locpad);
%
% vols{N,K} = load_md(md, nifti_modalities); if nargout == 1
% fnvols{N,K} = load_md(md, nifti_modalities, func); 
% [patches, layeridx, volidx, niifiles] = load_md(md, nifti_modalities, patchSize, locations, ...); if nargout <= 4
% [patches, layeridx, volidx, niifiles, vols] = load_md(md, nifti_modalities, patchSize, locations, ...); if nargout == 5

    narginchk(2, 7);
    nargoutchk(1, 5);
    md = varargin{1};
    % assuming matfile that has "md" in it, load it if necessary
    if ischar(md) 
        load(md);
    end
    modalities = varargin{2};
    
    isc = iscell(modalities);
    ismatfile = false;
    if ~isc
        f = md.getModality(modalities, 1);
        ismatfile = strcmp(f(end-3:end), '.mat');
    end
    
    if ismatfile % if matfile modality.

        % load matfile pointers
        mfptrs = arrayfunc(@(x) matfile(md.getModality(modalities, x)), 1:md.getNumSubjects);

        switch nargout 
            case 1
                varargout = {mfptrs};
            case 2
                vols = load_matfiles(mfptrs, varargin{3:end});
                varargout = {mfptrs, vols};
            case [3, 4]
                [patches, layeridx, volidx] = load_matfiles(mfptrs, varargin{3:end});
                varargout = {patches, layeridx, volidx, mfptrs};
                varargout = varargout(1:nargout);
            case 5
                [patches, layeridx, volidx, vols] = load_matfiles(mfptrs, varargin{3:end});
                varargout = {patches, layeridx, volidx, mfptrs, vols};
        end
         
    else % normal modalities.
        if ~iscell(modalities), modalities = {modalities}; end
        niifiles = cell(md.getNumSubjects(), numel(modalities));
        for v = 1:md.getNumSubjects()
            for m = 1:numel(modalities)
                niifiles{v, m} = md.getModality(modalities{m}, v);
            end
        end
        
        dopatches = numel(varargin) >= 4;
        
        if ~dopatches
            assert(nargout == 1);
            
            if nargin == 2
                varargout{1} = load_niifiles(niifiles);
            else
                assert(nargin == 3)
                varargout{1} = load_niifiles(niifiles, varargin{3});
            end
            
        elseif dopatches && nargout <= 4 % [patches, layeridx, volidx, niifiles]
            assert(numel(varargin) >= 3);
            [patches, layeridx, volidx] = load_niifiles(niifiles, varargin{3:end}); 
            varargout = {patches, layeridx, volidx, niifiles};
        else
            assert(dopatches && nargout == 5); % [patches, layeridx, volidx, niifiles, vols]
            [patches, layeridx, volidx, vols] = load_niifiles(niifiles, varargin{3:end}); 
            varargout = {patches, layeridx, volidx, niifiles, vols};
        end
    end
    
    varargout = varargout(1:nargout);
end
    
function varargout = load_matfiles(mfptrs, varargin)
% [vols{N,K}] = load_md(mfptrs, volumes2load); mfptrs can be matfile pointers of filenames
% [croppedvols{N,K}] = load_md(mfptrs, volumes2load, volcropsize); volcropsize is a cell so that cropvol = vol(volcropsize{:})
% [patches, layeridx, volidx] = load_md(mfptrs, volumes2load, patchSize, locations, <locpad>); 
% [patches, layeridx, volidx, vols] = load_md(mfptrs, volumes2load, patchSize, locations, <locpad>); 

    assert(ismember(nargin, [1, 2, 3, 4, 5]));
    if nargin == 4, varargin{4} = varargin{2} * 0; end
    nargoutchk(1, 4);
    if ~iscell(mfptrs), mfptrs = {mfptrs}; end
    if ischar(mfptrs{1}) % prepare matfiles if given strings
        mfptrs = cellfunc(@matfile, mfptrs);
    end
    
    dopatches = nargin >= 4;
    
    dovolload = ~dopatches || nargout == 4;
    if dovolload 
        cropvolsize = ifelse(nargin == 3, varargin(2), {});
        % cell of matfile pointers to volumes
        vols = subspacetools.matfiles2volumes(mfptrs, varargin{1}, cropvolsize{:});
    end
    
    if dopatches && dovolload % vols & patches
        [patches, layeridx, volidx] = load_volumes(vols, varargin{2:end});
        varargout = {patches, layeridx, volidx, vols};
        
    elseif ~dopatches && dovolload % just vols
        varargout = {vols};
        
    else
        assert(dopatches && ~dovolload) % just patches
        colcmd = @(v) subspacetools.matfile2patchcol(mfptrs, v, varargin{2:end});
        
        volnames = ifelse(iscell(varargin{1}), varargin{1}, varargin(1));
        patches = cell(numel(volnames), 1);
        for i = 1:numel(volnames)
            [patches{i}, layeridx, volidx] = colcmd(volnames{i});
        end
        
        patches = ifelse(iscell(varargin{1}), patches, patches{1});  % if volnames was a cell, keep patches as cell.
        varargout = {patches, layeridx, volidx};
    end
    
    varargout = varargout(1:nargout);
end

function varargout = load_volumes(volumes, varargin)
% vols = load_volumes(volumes{N,K})
% fnvols = load_volumes(volumes{N,K}, volfcn) % function to apply to volumes
% [patches, layeridx, volidx] = load_volumes(volumes, patchSize, locations, <locpad>) 
% [patches, layeridx, volidx, vols] = load_volumes(volumes, patchSize, locations, <locpad>) 
   
    narginchk(1, 4);
    if nargin == 3, varargin{3} = varargin{2} * 0; end
    nargoutchk(1, 4);
    vols = ifelse(iscell(volumes), volumes, {volumes});
    
    dopatches = nargin >= 3;
    
    if dopatches
        colcmd = @(vol) subspacetools.vols2patchcol(vol, varargin{:});
        patches = cell(1, size(vols, 2));
        for k = 1:size(vols, 2)
            [patches{k}, layeridx, volidx] = colcmd(vols(:, k));
        end 
        patches = ifelse(iscell(volumes), patches, patches{1}); % if volumes was a cell, keep patches as cell.
        varargout = {patches, layeridx, volidx};
        if nargout == 4
            varargout{4} = vols;
        end
        varargout = varargout(1:nargout);
    else
        if nargin == 1
            varargout = {vols};
        else
            varargout{1} = cellfun(varargin{2}, vols);
        end
    end
end
