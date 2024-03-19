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
bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);

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
fmin  = 5;
fmax = 180;
nfreq = 70;
frex = logspace(log10(fmin),log10(fmax),nfreq);

fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

TYPE_DB = 'default';
PATHS = struct();
PATHS.DataDir      = fullfile('W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results');
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


%%

% [E_speechlocked_all, E_speechlocked_avg_all, E_speechlocked_detr_all] = deal(cell(1,n_SUBJECTS));
%
% toDo = [1:4 6:n_SUBJECTS];
% for subj_i = toDo
%
%     % TF = [];
%     SUBJECT = SUBJECTS{subj_i};
%     fprintf('Starting computing TF results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
%     PATH_E = fullfile(PATHS.DataDir ,SUBJECT,'E.mat');
%     PATH_TF = fullfile(PATHS.DataDir ,SUBJECT,'TF.mat');
%
%     PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');
%
%     load(PATH_TF, 'E_speechlocked')
%     % save subjects
%     E_speechlocked_all{subj_i} = E_speechlocked;
%
% end

%% get events & MNI
% selected the three neurons, gets the highest channel coupled in
% a reigon of interest

toDo = 1 : n_SUBJECTS;% 1 : n_SUBJECTS;
[nSign_pairs_all, nPairs_all]  = deal(nan(1,n_SUBJECTS)) ;
Clusters_all = [];
PPC_mat_all = [];
PPCz_mat_all = [];
PPCmedianperm_mat_all = [];
ES_mat_all = [];
Phase_mat_all = [];
PairsLocation_MNI_all = [];
PLVTimeEvts_all = [];
for subj_i = toDo %
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)

    PATH_clusters = fullfile(PATHS.DataDir,SUBJECT,'ClustersPLV.mat');

    if isfile(PATH_clusters)
        fprintf(' Stacking clusters results in %s \n', PATH_clusters)
        load(PATH_clusters,'Clusters', 'nSign_pairs','n_pairs','PairsLocation_MNI', 'PPC_mat','PPCmedianperm_mat','PPCz_mat','ES_mat','PLVTimeEvts','Phase_mat');
        % put info about subjects id
        [Clusters.subj_id] = deal(SUBJECT);
        PairsLocation_MNI.subj_id = repmat(SUBJECT, height(PairsLocation_MNI),1);
        % stack Clusters
        Clusters_all = [Clusters_all Clusters(~isnan([Clusters.S_typeFRmod]))]; % eliminate nan firemod
        PairsLocation_MNI_all = [PairsLocation_MNI_all; PairsLocation_MNI];
        PPC_mat_all = [PPC_mat_all  PPC_mat];
        PPCz_mat_all = [PPCz_mat_all  PPCz_mat];
        PPCmedianperm_mat_all = [PPCmedianperm_mat_all  PPCmedianperm_mat];
        ES_mat_all = [ES_mat_all  ES_mat];
        Phase_mat_all = [Phase_mat_all Phase_mat];
        PLVTimeEvts_all = [PLVTimeEvts_all PLVTimeEvts];
        % grab information about significant pairs and n_pairs
        nSign_pairs_all(subj_i) = nSign_pairs;
        nPairs_all(subj_i) = n_pairs;

        fprintf(' Completed Stacking clusters results in %s \n', PATH_clusters)
        % clear variables for consistency
        clear Clusters nSign_pairs n_pairs PairsLocation_MNI PPC_mat PPCz_mat ES_mat PPCmedianperm_mat PLVTimeEvts Phase_mat
    else
        warning('Analysis still running: Clusters  is not available yet for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    end
end

nClusters = numel(Clusters_all);
nonan = ~isnan(nSign_pairs_all); % to manage the lack of data for patient 5



% Step 2.5. unwrap cell of cell arrray

Sign_Task = {};
%Sign_Task = [Sign_Task{:}];
PPC = [Clusters_all.PPC];
PPCz = [Clusters_all.PPCz];
ES = [Clusters_all.ES];

ClusterDur = [Clusters_all.TimeSpan];
ClusterCycles = [Clusters_all.nCycles];
ClusterZstat = [Clusters_all.Zstat];
ClusterFreqSpread =[Clusters_all.FreqSpread];
ClusternSign = [Clusters_all.nSign];
ClusterPhaseFlip=[Clusters_all.flip_polarity];

% this is a little bit more tricky!
ClusterPhase = [];
ClusterCentroid = [];
centerEvts = [];
ClusterOnOff = [];
ClusterE_MNI = [];
ClusterS_MNI = [];
ClusterPhaseinCluster = {Clusters_all.PhaseInCluster};
ClusterPhaseinCluster = [ClusterPhaseinCluster{:}];
ClusterSubjects = [];
ClusterEAtlasLabel = {};
ClusterEvtTimes = [];
ClusterEAtlasLabel_HCMMP = {};

Cluster_Schannel = {};
for clus_i = 1 : nClusters
    ClusterCentroid =[ClusterCentroid; Clusters_all(clus_i).Centroid];
    ClusterPhase =[ ClusterPhase; Clusters_all(clus_i).Phase];
    centerEvts = [centerEvts; Clusters_all(clus_i).centerEvts'];
    ClusterOnOff = [ClusterOnOff; [Clusters_all(clus_i).On' Clusters_all(clus_i).Off']];
    ClusterE_MNI = [ ClusterE_MNI; repmat([Clusters_all(clus_i).E_MNI_X,Clusters_all(clus_i).E_MNI_Y, Clusters_all(clus_i).E_MNI_Z],ClusternSign(clus_i),1)];
    ClusterS_MNI = [ ClusterS_MNI; repmat([Clusters_all(clus_i).S_MNI_X,Clusters_all(clus_i).S_MNI_Y, Clusters_all(clus_i).S_MNI_Z],ClusternSign(clus_i),1)];
    ClusterSubjects = [ClusterSubjects; repmat(Clusters_all(clus_i).subj_id,ClusternSign(clus_i),1)];
    Sign_Task = [Sign_Task ;Clusters_all(clus_i).window_task(:)];
    ClusterEAtlasLabel = [ClusterEAtlasLabel; repmat(Clusters_all(clus_i).E_atlas_label_Destrieux,ClusternSign(clus_i),1)];
    ClusterEvtTimes = [ClusterEvtTimes; repmat(median(Clusters_all(clus_i).TimeEvts.Evts),ClusternSign(clus_i),1)];
    ClusterEAtlasLabel_HCMMP = [ClusterEAtlasLabel_HCMMP; repmat(Clusters_all(clus_i).E_HCPMMP1_label_1,ClusternSign(clus_i),1)];
    Cluster_Schannel = [Cluster_Schannel; repmat(Clusters_all(clus_i).S_channel,ClusternSign(clus_i),1)];

end
ClusterSubjects = cellstr(ClusterSubjects);
nClusters_all = size(ClusterOnOff,1);

ClusterCentroidBand = nan(1,nClusters_all);
for clus_i = 1 : nClusters_all
    ClusterCentroidBand(clus_i) = find(ClusterCentroid(clus_i,2) < bands.fends,1);
end

meanEvts= median(ClusterEvtTimes);
clust_res = 0.005;
T = -2.5 : clust_res : 2;
%T = (min(ClusterEvtTimes(:))-.1) :clust_res: (max(ClusterEvtTimes(:)) + .1);
cluster_prob = zeros(nClusters_all,numel(T));
for ii = 1 : nClusters_all
    %     cluster_prob(ii,T<= (ClusterOnOff(ii,2) - meanSpeechOnset) & T >= (ClusterOnOff(ii,1) - meanSpeechOnset)) = 1;
    cluster_prob(ii,T<= (ClusterOnOff(ii,2)) & T >= (ClusterOnOff(ii,1))) = 1;
end

Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};

SignTask_id = zeros(5,nClusters_all);
for ii = 1 : nClusters_all
    tmp = Sign_Task{ii};
    for kk = 1:numel(tmp)
        SignTask_id(find(strcmpi(Task_labels,tmp{kk})),ii) = 1;
    end
end

bootsp_chance = 5*[0.34168 0.296690246157020]/100; % estimated before;
bootsp_chancesign = bootsp_chance(1) + bootsp_chance(2)*randn(1, 1000);
if all((bootsp_chancesign*sum(nPairs_all,[],'omitnan')) <= sum(nSign_pairs_all,[],'omitnan')  )
    pval_chancesign = 1/1000;
else
    pval_chancesign = sum(sum(nSign_pairs_all,[],'omitnan') >= (bootsp_chancesign*sum(nPairs_all,[],'omitnan')))/(1000 + 1);
end

%%

% plot

high_gamma_idx = frex >= bands.fstarts(6) &  frex <= bands.fends(6);
beta_idx = frex >= bands.fstarts(3) & frex <= bands.fends(4);
load('W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new\default\group-level\MNI_group.mat');
load('W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new\default\group-level\EVENTS_TASK_group.mat');
PATH_FIGURES = fullfile(PATH_OUTPUT,'Figures/TFwav');

for subj_i = 6% [21:25 16:20]
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Opening TF wavelet plots for SUBJECT %s, %d-%d \n', SUBJECT,subj_i,n_SUBJECTS)
    PATH_TF_CUE = fullfile(PATHS.DataDir ,SUBJECT,'TF_cue.mat');
    load(PATH_TF_CUE, 'E_cuelocked')
    PATH_TF_SPEECH = fullfile(PATHS.DataDir ,SUBJECT,'TF_speech.mat');
    load(PATH_TF_SPEECH, 'E_speechlocked')
    %
    evts = evts_all(subj_i);
    MNI_E = MNI_all(subj_i);

    %STG_list  = find(contains(MNI_E.E.atlas_label_Destrieux,"G_temp_sup"));


    % cue locked, speech locked and high gamma both curves

    for figx = 1 : height(MNI_E.E)
        figure('renderer','painters','position',[300 100 1200 300*numel(E_cuelocked)],'visible','off')
        tiledlayout(numel(E_cuelocked),3)

        if  true %ccccc~ismember(subj_i,[21:25 16:18])
            for session_id = 1 : numel(E_cuelocked)

                nexttile

                % plot cue-locked
                TF = E_cuelocked{1,session_id}.TF;
                time_cue = E_cuelocked{1,session_id}.time{1,1};
                baseline_time =  - [.6 .1];
                TFbase = squeeze(mean(TF(:,:,time_cue <= baseline_time(2) & time_cue >= baseline_time(1)),3,'omitnan'));
                TF_basenorm = 20*log10(TF./TFbase);
                TF_cue = squeeze(mean(TF_basenorm(figx,:,:),1,'omitnan'));

                %             h = pcolor(time_cue,frex,TF_cue);
                %             h.EdgeColor = 'none';
                imagesc(time_cue, frex, TF_cue)
                hold on
                plot((median(evts.cueon_evts) - median(evts.cueon_evts))*ones(1,2), [0 160], '--r')
                plot((median(evts.cueoff_evts) - median(evts.cueon_evts))*ones(1,2), [0 160], '--r')
                plot((median(evts.speechon_vts) - median(evts.cueon_evts))*ones(1,2), [0 160], '--r')
                plot((median(evts.speechoff_vts) - median(evts.cueon_evts))*ones(1,2), [0 160], '--r')
                xlabel(' Time [s] ')
                ylabel(' Frequency [Hz] ')
                set(gca,'YDir','normal')
                colormap(linspecer)
                colorbar
                clim([-1.5 1.5])
                box off
                ylim([5 160])
                xlim([-1 2])
                set(gca,'YScale','log')


                % plot speech-locked
                nexttile
                TF = E_speechlocked{1,session_id}.TF;
                time_speech = E_speechlocked{1,session_id}.time{1,1};
                baseline_time =  median(evts.cueon_evts) - [.6 .1];
                TFbase = squeeze(mean(TF(:,:,time_speech <= baseline_time(2) & time_speech >= baseline_time(1)),3,'omitnan'));
                TF_basenorm = 20*log10(TF./TFbase);
                TF_speech = squeeze(mean(TF_basenorm(figx,:,:),1,'omitnan'));

                %             h = pcolor(time_speech,frex,TF_speech);
                %             h.EdgeColor = 'none';
                imagesc(time_speech, frex, TF_speech)
                hold on
                plot(median(evts.speechon_vts)*ones(1,2), [0 160], '--r')
                plot(median(evts.speechoff_vts)*ones(1,2), [0 160], '--r')
                plot(median(evts.cueon_evts)*ones(1,2), [0 160], '--r')
                plot(median(evts.cueoff_evts)*ones(1,2), [0 160], '--r')

                xlabel(' Time [s] ')
                ylabel(' Frequency [Hz] ')
                set(gca,'YDir','normal')
                colormap(linspecer)
                colorbar
                clim([-1.5 1.5])
                box off
                ylim([5 160])
                xlim([-1 2])
                set(gca,'YScale','log')
                title(sprintf('Session %s', session_id))

                % plot curves
                nexttile
                high_gamma_cue = mean(TF_cue(high_gamma_idx,:),'omitnan');
                high_gamma_speech = mean(TF_speech(high_gamma_idx,:),'omitnan');
                plot(time_cue,high_gamma_cue)
                hold on
                plot(time_speech,high_gamma_speech)
                xlabel(' time [s] ')
                ylabel(sprintf('%s amplitude [dB]', char(bands.symbol(6))))
                box off
                xlim([-1 2])
                %ylim([-1 1])

                sgtitle(sprintf('Subject %s, Channel %s: %s',SUBJECT,char(MNI_E.E.electrode(figx)),char(MNI_E.E.atlas_label_Destrieux(figx))),'interpreter','none')
            end
            fprintf(' Saving electrode ERSP %d-%d \n', figx, height(MNI_E.E))
            figname = fullfile(PATH_FIGURES,sprintf('Subject-%s_Channel-%s_roiDestrieux-%s',SUBJECT,char(MNI_E.E.electrode(figx)),char(MNI_E.E.atlas_label_Destrieux(figx))));
            figname = strrep(figname,' ','');
            saveas(gcf,[figname,'.png'])
            saveas(gcf,[figname,'.fig'])
            saveas(gcf,[figname,'.svg'])
        end

        figure('renderer','painters','position',[300 100 800 300*numel(E_cuelocked)],'visible','off')
        tiledlayout(numel(E_cuelocked),2)

        for session_id = 1 : numel(E_cuelocked)

            TF = E_cuelocked{1,session_id}.TF;
            time_cue = E_cuelocked{1,session_id}.time{1,1};
            baseline_time =  - [.6 .1];
            TFbase = squeeze(mean(TF(:,:,time_cue <= baseline_time(2) & time_cue >= baseline_time(1)),3,'omitnan'));
            TF_basenorm = 20*log10(TF./TFbase);
            TF_cue = squeeze(mean(TF_basenorm(figx,:,:),1,'omitnan'));

            TF = E_speechlocked{1,session_id}.TF;
            time_speech = E_speechlocked{1,session_id}.time{1,1};
            baseline_time =  median(evts.cueon_evts) - [.6 .1];
            TFbase = squeeze(mean(TF(:,:,time_speech <= baseline_time(2) & time_speech >= baseline_time(1)),3,'omitnan'));
            TF_basenorm = 20*log10(TF./TFbase);
            TF_speech = squeeze(mean(TF_basenorm(figx,:,:),1,'omitnan'));

            beta_cue = mean(TF_cue(beta_idx,:),'omitnan');
            beta_speech = mean(TF_speech(beta_idx,:),'omitnan');
            high_gamma_cue = mean(TF_cue(high_gamma_idx,:),'omitnan');
            high_gamma_speech = mean(TF_speech(high_gamma_idx,:),'omitnan');

            nexttile

            plot(time_cue,beta_cue)
            hold on
            plot(time_speech,beta_speech)
            xlabel(' time [s] ')
            ylabel(sprintf('%s amplitude [dB]', char(bands.symbol(4))))
            box off
            xlim([-1 2])
            ylim([-1.5 1.5])
            nexttile


            plot(time_cue,high_gamma_cue)
            hold on
            plot(time_speech,high_gamma_speech)
            xlabel(' time [s] ')
            ylabel(sprintf('%s amplitude [dB]', char(bands.symbol(6))))
            box off
            xlim([-1 2])
            ylim([-1.5 1.5])



        end
        fprintf(' Saving electrode beta-gamma %d-%d \n', figx, height(MNI_E.E))
        figname = fullfile(PATH_FIGURES,sprintf('Subject-%s_Channel-%s_roiDestrieux-%s_freq-bands',SUBJECT,char(MNI_E.E.electrode(figx)),char(MNI_E.E.atlas_label_Destrieux(figx))));
        figname = strrep(figname,' ','');
        saveas(gcf,[figname,'.png'])
        saveas(gcf,[figname,'.fig'])
        saveas(gcf,[figname,'.svg'])
    end
    fprintf('Closing TF wavelet plots for SUBJECT %s, %d-%d \n', SUBJECT,subj_i,n_SUBJECTS)

end

