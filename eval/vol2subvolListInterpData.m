function vol2subvolListInterpData(interpMatFile, subjVolFile, subjMaskFile, atlSubvolLocFile, atlSubvolSize, outmatfile)
% from volume/subject data, for a given subvolume location and size in the
% atlas space, extract the appropriate original subvolume in subject space, as well as the forward
% and inverse R matrices.
% 
% optional (but used in patchSynthesis project: volIdxFile, a file with a variable "volIdx" that
% tells us *which* volume nrs to load in. This is because when preping the huge file matrix, we need
% to skip some volumes that don't have that appropriate modality.
%
% see also:md2subvolInterpData

    % input parsing
    narginchk(6, 6);
    if ischar(atlSubvolSize), atlSubvolSize = str2num(atlSubvolSize); end
    tic;
    
    % load file
    fid = fopen(atlSubvolLocFile);
    C = textscan(fid,'%d [%d,%d,%d]\n');
    fclose(fid);   
    inds = cat(1, C{1});
    atlSubvolLocs = cat(2, C{2:4});
    nLocs = numel(inds);

    % prepare output path
    mkdir(fileparts(outmatfile));
    
    % load in necessary variables
    interpMatData = load(interpMatFile, 'atlLoc2SubjSpace', 'atl2subjR', 'subj2atlR'); % cannot do partial loading :(
    subjVol = nii2vol(subjVolFile);
    subjMask = nii2vol(subjMaskFile);
    
    
    vi = verboseIter(1:nLocs, 2);
    while vi.hasNext()
        i = vi.next();
        [subjSubvol, subjMaskSubvol, R, G] = vol2subvolInterpData(interpMatData, subjVol, subjMask, atlSubvolLocs(i,:), atlSubvolSize);
        
        subjVarname = sprintf('subjSubvol_%d', inds(i));
        subjMaskVarname = sprintf('subjSubvolMask_%d', inds(i));
        rVarname = sprintf('R_%d', inds(i));
        gVarname = sprintf('Gamma_%d', inds(i));
        eval(sprintf('%s = subjSubvol;', subjVarname));
        eval(sprintf('%s = subjMaskSubvol;', subjMaskVarname));
        eval(sprintf('%s = R;', rVarname));
        eval(sprintf('%s = G;', gVarname));
        
        if sys.isfile(outmatfile)
            save(outmatfile, subjVarname, subjMaskVarname, rVarname, gVarname, '-append');
        else
            save(outmatfile, subjVarname, subjMaskVarname, rVarname, gVarname, '-v7.3');
        end
        
        % need to clear variables
        clear(subjVarname);
        clear(subjMaskVarname);
        clear(rVarname);
        clear(gVarname);
    end
    vi.close();
    