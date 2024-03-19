%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
PATH_OUTPUT = fullfile(PATH_PROJECT,'Supplementary_Analysis');

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

% PATH_FIGURES = ;

fmin  = 2;
fmax = 140;
nfreq = 60;

fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

TYPE_DB = 'default';
PATHS = struct();
PATHS.DataDir      = fullfile('W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results');
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


%%
evts_all = struct();
fs = 1000; % please check anytime

toDo = [1:4 6:n_SUBJECTS];

winHalfWidths_all = cell(1,n_SUBJECTS);
nSpikesWindow_all = cell(1,n_SUBJECTS);
for subj_i = toDo

    % TF = [];
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Starting computing PLV results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    PATH_PLV = fullfile(PATHS.DataDir ,SUBJECT,'PLV.mat');

    PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');
    % load PLV
    load(PATH_PLV)
    
    n_pairs = numel(PLV);
    winHalfWidths = nan(1,n_pairs);
    nSpikesWindow = nan(1,n_pairs);
    for pair_i = 1 : n_pairs 
            winHalfWidths(pair_i) = median(PLV(pair_i).store.all.winHalfWidths);
            nSpikesWindow(pair_i) = PLV(pair_i).store.nSpikes;       
    end

    winHalfWidths_all{subj_i} = winHalfWidths;
    nSpikesWindow_all{subj_i} = nSpikesWindow;
end

%% unwrap results
winHalfWidths_all_unwrap = cat(2,winHalfWidths_all{:});
nSpikesWindow_all_unwrap = cat(2,nSpikesWindow_all{:});

%% plot window properties

figure('Renderer','painters','Position',[300 300 300 600])
tiledlayout(2,1)

nexttile
h1 = raincloud_plot(winHalfWidths_all_unwrap'/fs + 0., 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
      'density_type','ks');
xlabel(' Half-width duration [s]')
%title(['Figure M6' newline 'Raincloud Plot: Reduced Smoothing, Kernel Density'])
%set(gca, 'XLim', [1 2]);
box off
nexttile
h1 = raincloud_plot(nSpikesWindow_all_unwrap', 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
      'density_type','ks');
%title(['Figure M6' newline 'Raincloud Plot: Reduced Smoothing, Kernel Density'])
%set(gca, 'XLim', [1 2]);
box off
xlabel(' # spikes PLV window ')

%%
figname = fullfile('Figures','windows_duration');
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'fig')