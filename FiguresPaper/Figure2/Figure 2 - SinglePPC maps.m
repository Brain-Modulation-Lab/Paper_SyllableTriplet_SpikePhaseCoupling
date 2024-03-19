clc
clear all
close all

%% Step 0: settings of the script
Figure_2_settings;

PATH_FIGURES = fullfile(PATH_OUTPUT,'Figures');

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
%cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-2.5 2.5];
cfg_fh.CMap = linspecer;
cfg_fh.Visible = 'off';
cfg_fh.smooth = true;

DB_TYPE = 'default';
fprintf('Printing pairs in all subjects')

% for loop across aptients
for subj_i = 6 : n_SUBJECTS
    fprintf('Printing pairs in subject-%s  (%d - %d) \n', SUBJECTS{subj_i}, subj_i,n_SUBJECTS)
    PATH_FIGURES_PATIENT = fullfile(PATH_FIGURES,DB_TYPE,SUBJECTS{subj_i});

    if ~isfolder(PATH_FIGURES_PATIENT)
        mkdir(PATH_FIGURES_PATIENT);
    end
    %% Step 1. gather 1 subject
    DB = get_SPCdata(SUBJECTS,subj_i, PATHS);
   

    % Step 2. Get single pairs with clusters
    if ~isempty(DB.Pairs.nPairs) | (DB.Pairs.nPairs == 0)
        SPCpair = get_SPCpair(DB, cfg);

        % Step 3. Plot single pair
        plotter_SPCpair(SPCpair, cfg_fh,PATH_FIGURES_PATIENT)
    end
    clear DB
end
fprintf('DONE! \n')


%%
% T = [-2.5 2];
% nPerms = 500;
% cfg.Baseline = [-Inf EvtTimes.Speech(2)];
% 
% SPCpair = get_SPCpairs(DB.Pairs,'PPCz', T, ...
%     Tres,cfg,'Speech','all','default',nPerms);
% 
% cfg_fh = struct();
% cfg_fh.Position = [300 300 700 500];
% cfg_fh.EvtTimes = EvtTimes.Speech;
% cfg_fh.EvtTypes = cfg.plot.EvtTypes;
% cfg_fh.CLim = [-.02 .18];
% cfg_fh.CMap = linspecer;
% 
% fh = plotter_SPCmap(SPCmap, cfg_fh)







