%% prepare the adni dataset

%% setup
% paths
atlasesDataSetup; % atlases
ADNI_PATH_ORIG = fullfile(GENERAL_DATA_PATH, '/ADNI/baselines/orig_brain');  % original files
ADNI_PATH_PROC = fullfile(SYNTHESIS_DATA_PATH, 'adni/proc'); % adni processing

% processing parameters
moveOrig2Proc = false; % move initial files from original folder to processing folder
getmdmethod = 'load'; % 'load', 'build', or 'none'
name = 'adni';
domdproc = false; % already done.
dsAmounts = 1:7;
intensityNorm = 255;
procdsRate = 5;

%% process
% if necessary, move files from original to srcpath from origpath
if moveOrig2Proc
    medicalDataset.files2folders(ADNI_PATH_ORIG, '\d*', ADNI_PATH_PROC);
end

% create a restorationmd
switch getmdmethod
    case 'build'
        md = restorationmd(dsAmounts, ADNI_PATH_PROC, SYNTHESIS_DATA_PATH, 'adni');
        
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
