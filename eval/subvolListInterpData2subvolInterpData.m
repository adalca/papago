function subvolListInterpData2subvolInterpData(path, sidsfile, subvolInd, outmatfile)
% path has a %s where the subjectid goes.
% sidsfile is a matfile with the 'sids' variable

    
    % load subject ids
    narginchk(4, 4);
    load(sidsfile, 'sids');
    if ischar(subvolInd), subvolInd = str2num(subvolInd); end
    
    % prepare variable names
    subjVarname = sprintf('subjSubvol_%d', subvolInd);
    rVarname = sprintf('R_%d', subvolInd);
    gVarname = sprintf('Gamma_%d', subvolInd);
    
    % prepare data
    nValidSubjects = numel(sids);
    subVols = cell(1, nValidSubjects);
    atl2SubjInterpMat = cell(1, nValidSubjects); % R in thesis
    subj2AtlInterpMat = cell(1, nValidSubjects); % Gamma in thesis
    
    % go through each subject
    for i = 1:nValidSubjects
        sid = sids{i};
        fullpath = sprintf(path, sid);
        
        % load necessary data for this subvol
        q = load(fullpath, subjVarname, rVarname, gVarname);
        subVols{i} = q.(subjVarname);
        atl2SubjInterpMat{i} = q.(rVarname);
        subj2AtlInterpMat{i} = q.(gVarname);
    end
    
    % save matfile
    tic;
    mkdir(fileparts(outmatfile));
    save(outmatfile, 'subVols', 'atl2SubjInterpMat', 'subj2AtlInterpMat', '-v7.3');
    toc;
