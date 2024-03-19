%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
PATH_OUTPUT = fullfile(PATH_PROJECT,'groupanalyses','Output');

PATH_UTILITIES =  fullfile(PATH_PROJECT,'Utilities');
addpath(genpath(PATH_UTILITIES))

prettify_plots;
ft_defaults;
bml_defaults;

% laod subjects
SUBJECTS = readtable(fullfile(PATH_PROJECT,'Subjects_list.txt'));
SUBJECTS = SUBJECTS.Subjects;

n_SUBJECTS = numel(SUBJECTS);

% parameters

cfg = set_configs("default-onlybeta");

%
% % load FR stats
% load(fullfile(PATH_RESULTS,"group-level","FRSTATS_group.mat"))
% FR_STATS = struct2table(FR_STATS);


% load cortex information
% load(fullfile(PATH_RESOURCES,'cortex_MNI.mat'),'BS1');


%
% FOI = [6 180];
% nfrex = 45;
% CYCLES = 5;
% TF_BUFFER = [0.5 0.5]; % add 3 sec before and after for wavelet decomposition
% TF_BASELINE_DUR = .75;
% TF_EVENT = "syl1_onset";

TYPE_DB = 'default';
% PATH_FIGURES = ;



fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

PATHS = struct();
PATHS.DataDir      = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures_new',TYPE_DB);


%% check DBS3029

for subj_i = 22
    
    TF = [];
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Starting computing PLV results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    PATH_E = fullfile(PATHS.DataDir ,SUBJECT,'E.mat');
    PATH_S = fullfile(PATHS.DataDir ,SUBJECT,'S.mat');
    
%     PATH_MNI = fullfile(PATH_RESULTS,SUBJECT,'MNI.mat');
    %PATH_STA = fullfile(PATH_RESULTS,SUBJECT,'STA.mat');
    
    PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');
    

    % load E
    load(PATH_E)
    
    % load S
    load(PATH_S)
    
    % compute PLV;
    PLV = check_E(E, S, cfg);
   
    

%     
%         STS.sel_unit = sel_units;
%         STS.n_units = numel(STS.label);
%         STS.n_channels = numel(E.label); 
    
    % save plv value
    disp("save plv results... ")
    mybeep
    %save(fullfile(PATHS.saveMatFiles,SUBJECT,'PLV_onlybeta.mat'), 'PLV',"-v7.3")

end

tEnd = toc(tStart);
disp(" Script ends ...")
fprintf(" Time elapsed %.2f \n",tEnd)
mybeep % play sound of end script
