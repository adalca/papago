function md2subvolInterpData(mdpath, interpMatMod, subjVolMod, atlSubvolLoc, atlSubvolSize, outmatfile, sidsFile)
% from a dataset stored as a medicalDataset object, for a given subvolume location and size in the
% atlas space, extract the appropriate original subvolume in subject space, as well as the forward
% and inverse R matrices.
% 
% optional (but used in patchSynthesis project: volIdxFile, a file with a variable "volIdx" that
% tells us *which* volume nrs to load in. This is because when preping the huge file matrix, we need
% to skip some volumes that don't have that appropriate modality.

    % input parsing
    narginchk(6, 7);
    if ischar(atlSubvolLoc), atlSubvolLoc = str2num(atlSubvolLoc); end
    if ischar(atlSubvolSize), atlSubvolSize = str2num(atlSubvolSize); end
    tic;
    
    % load md file 
    load(mdpath);
    
    % get subject nrs
    if exist('sidsFile', 'var')
        load(sidsFile, 'sids');
        [~, ~, volIdx] = intersect(sids, md.sids);
        assert(numel(volIdx) == numel(sids));
    else
        sids = md.sids;
        volIdx = 1:numel(sids);
    end
    nValidSubjects = numel(sids);
    
    % prepare data
    subVols = cell(1, nValidSubjects);
    atl2SubjInterpMat = cell(1, nValidSubjects); % R in thesis
    subj2AtlInterpMat = cell(1, nValidSubjects); % Gamma in thesis
    
    % load and parse data
    vi = verboseIter(volIdx, 2);
    while vi.hasNext()
        [vii, idx] = vi.next();
        
        % load in necessary variables
        interpMatFile = md.getModality(interpMatMod, vii);
        load(interpMatFile, 'atlLoc2SubjSpace', 'atl2subjR', 'subj2atlR'); % cannot do partial loading :(
        subjVol = md.loadVolume(subjVolMod, vii);
        
        % prep atlas size linear indexes. 
        % should be done outside the loop, but not sure of atlas size...
        if ~exist('atlsize', 'var')
            atlsize = size(atlLoc2SubjSpace{1});
            atlRange = arrayfunc(@(m,s) m:m+s-1, atlSubvolLoc, atlSubvolSize);
            atlRange = ndgrid2cell(atlRange{:});
            atlIdx = sub2ind(atlsize, atlRange{:});
        end
        
        % get subject subspace
        subjSubvolRange = srcVol2tgtPatch(atlSubvolLoc, atlSubvolSize, atlLoc2SubjSpace);
        subjSubvolRangeNd = ndgrid2cell(subjSubvolRange{:});
        
        % get the linear indexes for both atl and subject space, so as to extract the appropriate
        % part of atl2subjR and subj2atlR (extracting both is extremly large)
        %   note on orientation: atl2subjR is #SubjVox-by-#AtlVox, so it takes the atlas to the
        %   subject. In thesis notation, this is "R". subj2atlR is then Gamma.
        subjIdx = sub2ind(size(subjVol), subjSubvolRangeNd{:});
        R = atl2subjR(subjIdx(:), atlIdx(:));
        G = subj2atlR(atlIdx(:), subjIdx(:));
        
        % add results to cell
        subVols{idx} = subjVol(subjSubvolRange{:});
        atl2SubjInterpMat{idx} = R;
        subj2AtlInterpMat{idx} = G;
    end
    vi.close();
    toc;
        
    nfo = struct();
    nfo.subvolLoc = atlSubvolLoc;
    nfo.subvolSize = atlSubvolSize;
    
    % save matfile
    tic;
    mkdir(fileparts(outmatfile));
    save(outmatfile, 'subVols', 'atl2SubjInterpMat', 'subj2AtlInterpMat', 'nfo', '-v7.3');
    toc;
    