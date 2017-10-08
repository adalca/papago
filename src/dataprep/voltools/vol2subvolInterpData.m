function [subjSubvol, subjMaskSubvol, R, G, rangeMins] = ...
    vol2subvolInterpData(interpMatData, subjVol, subjMask, atlSubvolLoc, atlSubvolSize, outmatfile)
% from volume/subject data, for a given subvolume location and size in the
% atlas space, extract the appropriate original subvolume in subject space, as well as the forward
% and inverse R matrices.
% 
% optional (but used in patchSynthesis project: volIdxFile, a file with a variable "volIdx" that
% tells us *which* volume nrs to load in. This is because when preping the huge file matrix, we need
% to skip some volumes that don't have that appropriate modality.
%
% see also:md2subvolInterpData
%
% TODO: replace some of this vode with vol2subvolInterpmat() !!!

    % input parsing
    narginchk(5, 6);
    if ischar(atlSubvolLoc), atlSubvolLoc = str2num(atlSubvolLoc); end
    if ischar(atlSubvolSize), atlSubvolSize = str2num(atlSubvolSize); end
    tic;
    
    if ischar(subjVol) && sys.isfile(subjVol)
        subjVol = nii2vol(subjVol);
    end
    
    if ischar(subjMask) && sys.isfile(subjMask)
        subjMask = nii2vol(subjMask);
    end
    
    if ischar(interpMatData) && sys.isfile(interpMatData)
        interpMatData = load(interpMatData);
    end
        
    % prep atlas size linear indexes. 
    atlsize = size(interpMatData.atlLoc2SubjSpace{1});
    atlRange = arrayfunc(@(m,s) m:m+s-1, atlSubvolLoc, atlSubvolSize);
    atlRange = ndgrid2cell(atlRange{:});
    atlIdx = sub2ind(atlsize, atlRange{:});

    % get subject subspace
    [subjSubvolRange, rangeMins] = srcVol2tgtPatch(atlSubvolLoc, atlSubvolSize, interpMatData.atlLoc2SubjSpace);
    subjSubvolRangeNd = ndgrid2cell(subjSubvolRange{:});

    % get the linear indexes for both atl and subject space, so as to extract the appropriate
    % part of atl2subjR and subj2atlR (extracting both is extremly large)
    %   note on orientation: atl2subjR is #SubjVox-by-#AtlVox, so it takes the atlas to the
    %   subject. In thesis notation, this is "R". subj2atlR is then Gamma.
    subjIdx = sub2ind(size(subjVol), subjSubvolRangeNd{:});
    R = interpMatData.atl2subjR(subjIdx(:), atlIdx(:));
    % G = interpMatData.subj2atlR(atlIdx(:), subjIdx(:));
    warning('G not being read. FIXME');
    G = nan;

    nfo = struct();
    nfo.subvolLoc = atlSubvolLoc;
    nfo.subvolSize = atlSubvolSize;
    
    subjSubvol = subjVol(subjSubvolRange{:});
    subjMaskSubvol = subjMask(subjSubvolRange{:});
    
    % save matfile
    if exist('outmatfile', 'var')
        tic;
        mkdir(fileparts(outmatfile));
        save(outmatfile, 'subjSubvol', 'subjMaskSubvol', 'R', 'G', 'nfo', '-v7.3');
        toc;
    end
    