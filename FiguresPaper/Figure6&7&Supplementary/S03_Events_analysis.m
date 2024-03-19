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

toDo = [1:4 6:n_SUBJECTS];
for subj_i = toDo

    % TF = [];
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Starting computing PLV results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    PATH_PLV = fullfile(PATHS.DataDir ,SUBJECT,'ClustersPLV.mat');

    PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');
    % load PLV
    load(PATH_PLV, 'PLVTimeEvts')
    
    PLV_Evts = PLVTimeEvts.Evts;

    evts_all(subj_i).id = SUBJECT;

    % save information
    evts_all(subj_i).cue_duration = PLV_Evts(:,3) - PLV_Evts(:,2);
    evts_all(subj_i).speech_duration = PLV_Evts(:,5) - PLV_Evts(:,4);
    evts_all(subj_i).prespeech_duration = PLV_Evts(:,4) - PLV_Evts(:,3);
    evts_all(subj_i).cueon_evts = PLV_Evts(:,2);
    evts_all(subj_i).cueoff_evts = PLV_Evts(:,3);
    evts_all(subj_i).speechon_vts = PLV_Evts(:,4);
    evts_all(subj_i).speechoff_vts = PLV_Evts(:,5);
end

%% convert struct to array
cue_duration_unwrap = [];
speech_duration_unwrap = [];
prespeech_duration_unwrap = [];
iti_duration_unwrap = [];
cueon_evts_unwrap = [];
cueoff_evts_unwrap = [];
speechon_evts_unwrap = [];
speechoff_evts_unwrap = [];
for subj_i = 1 : n_SUBJECTS
    cue_duration_unwrap = [cue_duration_unwrap; evts_all(subj_i).cue_duration];
    speech_duration_unwrap = [speech_duration_unwrap; evts_all(subj_i).speech_duration];
    prespeech_duration_unwrap = [prespeech_duration_unwrap; evts_all(subj_i).prespeech_duration];
    cueon_evts_unwrap = [cueon_evts_unwrap; evts_all(subj_i).cueon_evts];
    cueoff_evts_unwrap = [cueoff_evts_unwrap; evts_all(subj_i).cueoff_evts];
    speechon_evts_unwrap = [speechon_evts_unwrap; evts_all(subj_i).speechon_vts];
    speechoff_evts_unwrap = [speechoff_evts_unwrap; evts_all(subj_i).speechoff_vts];
end
%%
% patient DBS3008 (number 6)
subj_i = 6;
data_toplot = [evts_all(subj_i).cueon_evts evts_all(subj_i).cueoff_evts evts_all(subj_i).speechon_vts evts_all(subj_i).speechoff_vts];
[~,idx_sort_speechdur] = sort(evts_all(subj_i).speech_duration,'ascend');
data_toplot_sort = data_toplot(idx_sort_speechdur,:);
nevents = size(data_toplot,1);
figure('renderer','painters')
plot([data_toplot_sort(:,1) data_toplot_sort(:,2)]',repmat(1:nevents,2,1),'linewidth',1.2,'color','k')
hold on
plot([data_toplot_sort(:,3) data_toplot_sort(:,4)]',repmat(1:nevents,2,1),'linewidth',1.2,'color','k')
box off
ylabel(' ordered trials id [#]')
xlabel(' Time [s] ')
xlim([-4 3])

% patient percentage cue or speech
T_trial = -4 : 0.05 : 3;
[trial_perc_cue, trial_perc_speech] = deal(nan(n_SUBJECTS,numel(T_trial)));

for subj_i = 1 : n_SUBJECTS
    if ~isempty(evts_all(subj_i).speechon_vts)
        trial_perc_speech(subj_i,:) = mean(T_trial >= evts_all(subj_i).speechon_vts & T_trial <= evts_all(subj_i).speechoff_vts);
        trial_perc_cue(subj_i,:) = mean(T_trial >= evts_all(subj_i).cueon_evts & T_trial <= evts_all(subj_i).cueoff_evts);
    end
end

speechdur = arrayfun(@(x) mean(x.speech_duration),evts_all);
[~,idx_sort_speechdur] = sort(speechdur,'descend');

color_orangemap = [linspace(1,1,45)', linspace(1,0.55,45)', linspace(1,0,45)'];
color_blackmap = [linspace(1,0,45)', linspace(1,0,45)', linspace(1,0,45)'];

figure

imagesc(T_trial, (1:n_SUBJECTS)-1, trial_perc_speech(idx_sort_speechdur,:))
colormap(color_orangemap);
xlabel(' Time [s] ')
ylabel(' ordered patients [#] ')
cb = colorbar;
ylabel(cb,'trials [%]')
box off
xlim([-4.5 3.5])
 ylim([1 24])
 yticks(1:5:25)
% yticks(2:25)
figname = fullfile('Figures','events_subj_speech_trials');
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'fig')
saveas(gcf,figname,'svg')

figure
imagesc(T_trial, (1:n_SUBJECTS)-1, trial_perc_cue(idx_sort_speechdur,:))
colormap(color_blackmap)
xlabel(' Time [s] ')
ylabel(' ordered patients [#] ')
ylabel(cb,'AP')
cb = colorbar;
ylabel(cb,'trials [%]')
box off
ylim([2 25])

xlim([-4.5 3.5])
 ylim([1 24])
 yticks(1:5:25)
 figname = fullfile('Figures','events_subj_cue_trials');
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'fig')
saveas(gcf,figname,'svg')

%% plot intervals

% plot duration
figure('Renderer','painters','Position',[300 300 300 900])
tiledlayout(3,1)

nexttile
h1 = raincloud_plot(cue_duration_unwrap', 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
      'density_type','ks');
xlabel(' cue duration [s]')
%title(['Figure M6' newline 'Raincloud Plot: Reduced Smoothing, Kernel Density'])
set(gca, 'XLim', [1 2]);
box off
nexttile
h1 = raincloud_plot(speech_duration_unwrap', 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
      'density_type','ks');
%title(['Figure M6' newline 'Raincloud Plot: Reduced Smoothing, Kernel Density'])
%set(gca, 'XLim', [1 2]);
box off
xlabel(' speech duration [s]')

nexttile
h1 = raincloud_plot(prespeech_duration_unwrap', 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
      'density_type','ks');
%title(['Figure M6' newline 'Raincloud Plot: Reduced Smoothing, Kernel Density'])
%set(gca, 'XLim', [1 2]);
box off
xlabel(' prespeech duration [s]')

figname = fullfile('Figures','events_duration');
saveas(gcf,figname,'png')
saveas(gcf,figname,'svg')
saveas(gcf,figname,'fig')




figure('Renderer','painters','Position',[300 300 600 300])
h1 = raincloud_plot(cueon_evts_unwrap, 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(cueoff_evts_unwrap, 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3 = raincloud_plot(speechoff_evts_unwrap, 'box_on', 1, 'color', [.6 .6 .6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
%legend([h1{1} h2{1} h3{1}], {'Group 1', 'Group 2'});
%title(['Figure M7' newline 'A) Dodge Options Example 1']);
%set(gca,'XLim', [0 40], 'YLim', [-.075 .15]);

box off
figname = fullfile('Figures','events_timing');
saveas(gcf,figname,'png')
saveas(gcf,figname,'svg')
saveas(gcf,figname,'fig')





%

tEnd = toc(tStart);
disp(" Script ends ...")
fprintf(" Time elapsed %.2f \n",tEnd)
mybeep % play sound of end script
