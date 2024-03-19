%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
PATH_PROJECT = 'Y:\Users\MV1019\PhaseLocking';
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
PATHS.DataDir      = fullfile('Y:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results');
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


PHONEME_MAPPING = readtable(fullfile(PATH_PROJECT,'phoneme_mapping.csv'),NumHeaderLines=1);
PHONEME_TARGET = PHONEME_MAPPING{:,1};
PHONEME_PROD = PHONEME_MAPPING{:,2};
% ignore the 10th row
PHONEME_TARGET(10,:) = [];
PHONEME_PROD(10,:) = [];

%
evts_all = struct();
allevts_speechdur = [];
toDo = [1:4 6:n_SUBJECTS];
nTrials = nan(n_SUBJECTS,1);
nCorrTrials = nan(n_SUBJECTS,1);
PercCorrTrials = nan(n_SUBJECTS,1);
PercCorrTargetTrials= nan(n_SUBJECTS,1);
for subj_i = toDo

    % TF = [];
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Starting computing PLV results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    PATH_PLV = fullfile(PATHS.DataDir ,SUBJECT,'ClustersPLV.mat');
    PATH_E = fullfile(PATHS.DataDir ,SUBJECT,'E.mat');
    PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');
    % load PLV
    load(PATH_PLV, 'PLVTimeEvts')

    % load accuracy info
    produced_phoneme = bml_annot_read_tsv(fullfile(PATH_ANNOT,[SUBJECT,'_produced_phoneme.txt']));
    produced_triplet = bml_annot_read_tsv(fullfile(PATH_ANNOT,[SUBJECT,'_produced_triplet.txt']));
    % syllid 
    trialid = unique(produced_phoneme.trial_id);
    sessid = unique(produced_phoneme.session_id);

    tmp = [];
    for sess_i = 1 : numel(sessid)
        for trial_i = 1 : numel(trialid)
            tmp_acc = produced_phoneme.accuracy(produced_phoneme.trial_id == trialid(trial_i) & produced_phoneme.session_id == sessid(sess_i));
            tmp = [tmp; all(strcmpi(tmp_acc,'accurate'))];

        end
    end
    nCorrTrials(subj_i) = sum(tmp);
    nTrials(subj_i) = numel(tmp);
    PercCorrTrials(subj_i) = mean(tmp)*100;


    % define wheter triplet is on target or not [phnetic_target_accuracy]

    

    list_triplets_id = produced_triplet.trial_id;
    list_triplets_session = produced_triplet.session_id;
    phonetic_ontarget_ = [];
    %phonetic_ontarget_perc_ = [];
    
    for tripl_i = 1 : height(produced_triplet)
        idx_phoneme = produced_phoneme.trial_id == list_triplets_id(tripl_i) &  produced_phoneme.session_id == list_triplets_session(tripl_i);
        phonemes_in_tripl = produced_phoneme(idx_phoneme,:);
        phonemes_in_tripl_target = phonemes_in_tripl{:,'stim'};
        phonemes_in_tripl_prod = phonemes_in_tripl{:,'phonetic_code'};
        phonemes_ontarget = cellfun(@(x,y) ismember(x, PHONEME_PROD(ismember(PHONEME_TARGET,y))), phonemes_in_tripl_prod,phonemes_in_tripl_target);
        produced_phoneme.phonetic_ontarget(idx_phoneme) = phonemes_ontarget;  

        is_ontarget =  all(phonemes_ontarget);
        perc_ontarget =  mean(phonemes_ontarget)*100;

        phonetic_ontarget_ = [phonetic_ontarget_;is_ontarget];
        produced_triplet.phonetic_ontarget(tripl_i) = is_ontarget;
        %phonetic_ontarget_perc_ = [phonetic_ontarget_perc_;perc_ontarget];
        produced_triplet.phonetic_ontarget_perc(tripl_i) = perc_ontarget;


    end
    

    PercCorrTargetTrials(subj_i) = mean(phonetic_ontarget_,'omitnan')*100;

    
    % PLV_Evts = [PLVTimeEvts.Speech.time];
    % 
    % evts_all(subj_i).id = SUBJECT;
    % 
    % % save information
    % evts_all(subj_i).cue_duration = PLV_Evts(:,3) - PLV_Evts(:,2);
    % evts_all(subj_i).speech_duration = PLV_Evts(:,5) - PLV_Evts(:,4);
    % evts_all(subj_i).cueon_evts = PLV_Evts(:,2);
    % evts_all(subj_i).cueoff_evts = PLV_Evts(:,3);
    % evts_all(subj_i).speechon_vts = PLV_Evts(:,4);
    % evts_all(subj_i).speechoff_vts = PLV_Evts(:,5);
    % 
    % allevts_speechdur = [allevts_speechdur; evts_all(subj_i).speech_duration ];


    % update anntoaiton table...
    bml_annot_write_tsv(produced_phoneme,fullfile(PATH_ANNOT,[SUBJECT,'_produced_phoneme.txt']));
    bml_annot_write_tsv(produced_triplet,fullfile(PATH_ANNOT,[SUBJECT,'_produced_triplet.txt']));
end
%%
meanevts = arrayfun(@(x) mean(x.speech_duration,'omitnan'), evts_all,'uni',true)';
stdevts = arrayfun(@(x) std(x.speech_duration,[],'omitnan'), evts_all,'uni',true)';

meanevts_pool = mean(allevts_speechdur,'omitnan');
stdevts_pool = std(allevts_speechdur,[],'omitnan');

%%

figure('Renderer','painters','Position',[300 300 200 300])
imagesc(PercCorrTargetTrials(~isnan(PercCorrTargetTrials)))
colormap(flipud(gray))
cb = colorbar;
ylabel(cb,'CVCVCV accuracy [%]')

%legend([h1{1} h2{1} h3{1}], {'Group 1', 'Group 2'});
%title(['Figure M7' newline 'A) Dodge Options Example 1']);
%set(gca,'XLim', [0 40], 'YLim', [-.075 .15]);


%
box off
figname = fullfile('Figures','accuracy_perf');
saveas(gcf,figname,'png')
saveas(gcf,figname,'svg')
saveas(gcf,figname,'fig')

