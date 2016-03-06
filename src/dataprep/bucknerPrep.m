%% prepare the buckner dataset

%% setup
% paths
atlasesDataSetup; % atlases
BUCKNER_PATH_ORIG = fullfile(BUCKNER_ATLAS_ORIG.PATH, 'orig'); % original files
BUCKNER_PATH_PROC = fullfile(SYNTHESIS_DATA_PATH, 'buckner/proc'); % buckner processing 

% processing parameters
moveOrig2Proc = false; % move initial files from original folder to processing folder
getmdmethod = 'load'; % 'load', 'build', or 'none'
name = 'buckner';
domdproc = false; % already done.
dsAmounts = 2:7;
intensityNorm = 255;
procdsRate = 5;

%% process
% if necessary, move files from original to srcpath from origpath
if moveOrig2Proc
    medicalDataset.files2folders(BUCKNER_PATH_ORIG, 'buckner\d*', BUCKNER_PATH_PROC);
end

% create buckner medicalDatasets for data.
switch getmdmethod
    case 'build'
        md = restorationmd(dsAmounts, BUCKNER_PATH_PROC, SYNTHESIS_DATA_PATH, name);
        
    case 'load'
        % get latest file
        fnames = fullfile(SYNTHESIS_DATA_PATH, name, 'md', [sys.usrname, '_restor_md_*']);
        md = loadmd(fnames);

    case 'none'
        clear md
    otherwise
        error('unknown md loading method');
end

% process md 
if domdproc
    processmd(md, procdsRate, intensityNorm, BUCKNER_ATLAS_MODS);
end
