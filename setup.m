%% Setup the subspace project
% This file loads important variables, paths and libraries for the subspace project. 
% It is meant to be run once at the beginning of working with the project.

%% Settings
warning off backtrace; % turn off backtrace for warnings.
warning off verbose;

%% PATH Settings
if ispc % adrian
    TOOLBOXES_PATH = 'C:/Users/adalca/Dropbox (Personal)/MATLAB/toolboxes';
    SYNTHESIS_DATA_PATH = 'D:/Dropbox (MIT)/Research/patchSynthesis/data/'; % synthetic processed data
    GENERAL_DATA_PATH = 'D:/Dropbox (MIT)/Research/generalData'; % original data
    OUTPUT_PATH = 'D:/Dropbox (MIT)/Research/patchSynthesis/output/subspace';
    
else
    [~, whoami] = system('whoami');
    spl = strsplit(whoami, '\');
    usrname = strtrim(spl{end}); 
    switch usrname
        case 'klbouman' % katie's local mac
            TOOLBOXES_PATH = '/Users/klbouman/Research/medicalInpainting/';
            SYNTHESIS_DATA_PATH = '/Users/klbouman/Dropbox (MIT)/patchSynthesis/data/';
            GENERAL_DATA_PATH = '/Users/klbouman/Dropbox (MIT)/';
            OUTPUT_PATH = '/Users/klbouman/Research/medicalInpainting/dump/';
        
        case 'adalca'
            TOOLBOXES_PATH = '/data/vision/polina/users/adalca/patchSynthesis/subspace/latest/toolboxes/';
            SYNTHESIS_DATA_PATH = '/data/vision/polina/scratch/adalca/patchSynthesis/data/';
            GENERAL_DATA_PATH = SYNTHESIS_DATA_PATH;
            OUTPUT_PATH = '/data/vision/polina/scratch/adalca/patchSynthesis/output/';
    end
        
end

%% add toolboxes and project code
addpath(genpath(TOOLBOXES_PATH));
addpath(genpath(pwd));

%% prepare a set of locations
LOC_VENTRICLE_EDGE = [56, 72, 92]; % for the original scale
LOC_LEFT_CORTEX = [40, 74 114];

%% set up data paths

% setup atlasses
atlasesDataSetup;

