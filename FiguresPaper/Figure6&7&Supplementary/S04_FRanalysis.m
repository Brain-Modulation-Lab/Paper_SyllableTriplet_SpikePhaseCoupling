%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
if ispc
    PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
else
    PATH_PROJECT = '/Volumes/Nexus4/Users/MV1019/PhaseLocking';
end
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
PATH_OUTPUT = fullfile(PATH_PROJECT,'Supplementary_Analysis');
PATH_DISTAL = 'Z:\Resources\DISTAL-Atlas';
PATH_NIFTI = fullfile(PATH_OUTPUT,'Nifti');

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
cfg = set_configs('default');

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
PATHS.saveMatFiles = fullfile(PATH_PROJECT,'groupanalyses','Output','Results_new',TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;


% build information abotu unti types used in the analysis
% PairsLocation_MNI_all = [];
% toDo = [1:4 6:n_SUBJECTS];
% for subj_i = toDo
%
%     % TF = [];
%     SUBJECT = SUBJECTS{subj_i};
%       fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
%     PATH_clusters = fullfile(PATHS.saveMatFiles ,SUBJECT,'ClustersPLV.mat');
%
%
%     if isfile(PATH_clusters)
%         fprintf(' Stacking clusters results in %s \n', PATH_clusters)
%         load(PATH_clusters,'PairsLocation_MNI');
%         % put info about subjects id
%         PairsLocation_MNI.subj_id = repmat(SUBJECT, height(PairsLocation_MNI),1);
%         % stack location
%         PairsLocation_MNI_all = [PairsLocation_MNI_all; PairsLocation_MNI];
%     end
% end

%%

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

    PATH_clusters = fullfile(PATHS.saveMatFiles,SUBJECT,'ClustersPLV.mat');

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



unit_types = [-1 0 1 2];
unit_types_labels= {'Neg','None','Excit','Mixed'};
n_unit_types_labels = [40, 69 85 25];

%%

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
ClusterS_FRMod = [];
ClusterSChannel = [];
ClusterEChannel = [];
ClusterSunitGrade = {};



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
    %ClusterEvtTimes = [ClusterEvtTimes; repmat(median(Clusters_all(clus_i).TimeEvts.Evts),ClusternSign(clus_i),1)];
    ClusterEAtlasLabel_HCMMP = [ClusterEAtlasLabel_HCMMP; repmat(Clusters_all(clus_i).E_HCPMMP1_label_1,ClusternSign(clus_i),1)];
    ClusterS_FRMod = [ClusterS_FRMod;  repmat(Clusters_all(clus_i).S_typeFRmod,ClusternSign(clus_i),1)];
    ClusterSChannel = [ClusterSChannel; repmat(string(Clusters_all(clus_i).S_channel),ClusternSign(clus_i),1)];
    ClusterEChannel = [ClusterEChannel; repmat(string(Clusters_all(clus_i).E_channel),ClusternSign(clus_i),1)];
    ClusterSunitGrade = [ClusterSunitGrade; repmat(string(Clusters_all(clus_i).S_unitGrade),ClusternSign(clus_i),1)];
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
%% give basic stats
PairLoc = [PairsLocation_MNI_all.S_MNI_X PairsLocation_MNI_all.S_MNI_Y PairsLocation_MNI_all.S_MNI_Z];
Unit_id = [PairsLocation_MNI_all.S_channel PairsLocation_MNI_all.E_channel PairsLocation_MNI_all.subj_id];
Unit_FRmod = PairsLocation_MNI_all.S_typeFRmod;
Unit_id = Unit_id(all(~isnan(PairLoc),2),:);
Unit_FRmod = Unit_FRmod(all(~isnan(PairLoc),2),:);

PairLoc = PairLoc(all(~isnan(PairLoc),2),:);

% number of pairs with at least one clusters
figure('renderer','painters','defaultAxesFontSize', 12, 'Position', [300 300 800 400])
bar(1,sum(nSign_pairs_all,[],'omitnan')/sum(nPairs_all,[],'omitnan')*100, "FaceColor", [.8 .8 .8] , 'FaceAlpha', 1, 'BarWidth',1);
hold on
bar((1 : sum(nonan)) + 4,nSign_pairs_all(nonan)./nPairs_all(nonan)*100, "FaceColor", [.8 .8 .8], 'FaceAlpha', .4);
xticks([1 (1 : sum(nonan)) + 4])
xticklabels([ 'All',SUBJECTS(nonan)'])
ylabel(" Sign. pairs [%] ")
box off
figname = fullfile(PATHS.saveFigures,'Sign_pairs_bar');
saveas(gcf,figname,'fig')
exportgraphics(gcf,[figname,'.png'])
exportgraphics(gcf,[figname,'.eps'])


% numbers of clusters per pair
nSigns = unique(ClusternSign);
figure('renderer','painters','position',[400 400 400 300])
hold on
arrayfun(@(x) bar(x,sum(ClusternSign == x),"FaceColor", [.8 .8 .8], 'FaceAlpha', .4),nSigns)
xlabel(" # clusters per pair ")
ylabel(" # pairs ")
xticks(nSigns)
% axes('Position',[.5 .5 .4 .4])
% box on
% arrayfun(@(x) bar(x,sum(ClusternSign > 1),"FaceColor", [.8 .8 .8], 'FaceAlpha', .4),nSigns)
text(4,1500, sprintf(['   %1.3f ' char(177)  ' %1.3f [a.u.] '], mean(ClusternSign) , std(ClusternSign)))

figname = fullfile(PATHS.saveFigures,'clusters_per_pair');

saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'fig')

% number of units with at least one clusters
PairsLocation_MNI_all.subj_id = string(PairsLocation_MNI_all.subj_id);
[uniqueUnits,ia_uniqueUnits,ix_uniqueUnits] = unique([ PairsLocation_MNI_all.subj_id PairsLocation_MNI_all.S_channel],'rows');
nuniqueUnits_perunittype = arrayfun(@(x) sum(PairsLocation_MNI_all.S_typeFRmod(ia_uniqueUnits) == x),  unit_types);
uniqueUnits_nPairs = accumarray(ix_uniqueUnits,1);
uniqueUnits_nClusters = arrayfun(@(x) sum(PairsLocation_MNI_all.nSign(ix_uniqueUnits == x)), unique(ix_uniqueUnits));
uniqueUnits_oneClust =  mean(uniqueUnits_nClusters >0)*100;
uniqueUnits_nClusters_bypair =  uniqueUnits_nClusters./uniqueUnits_nPairs;

figure('Renderer','painters','Position',[300 300 800 300])
tiledlayout(1,2)
nexttile()
raincloud_plot(uniqueUnits_nClusters_bypair,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(uniqueUnits_nClusters_bypair),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f '],median(uniqueUnits_nClusters_bypair),iqr(uniqueUnits_nClusters_bypair)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel('Clusters/pair per unit [a.u]')
ylabel(" PDF [a.u.] ")
yticks([0 2 4])
xlim([0 1])
nexttile()
raincloud_plot(uniqueUnits_nClusters,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(uniqueUnits_nClusters),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f '],median(uniqueUnits_nClusters),iqr(uniqueUnits_nClusters)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel('# Clusters per unit')
ylabel(" PDF [a.u.] ")
yticks([0 2 4])
xlim([0 100])
axes('Position',[.7 .7 .2 .2])
box on
raincloud_plot(uniqueUnits_nClusters,"box_on",1,"color",[.8 .8 .8],"alpha",.3,'sizedata',14)
xlim([0 8])
ylim([-0.05 -0.02])
figname = fullfile(PATHS.saveFigures,'clusters_per_pair');

saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'fig')

%% spatial location distribution of cell types (+, - , +-, x) (Create nifti file)
res_sphere = 1;

STN_nii = ea_load_nii(fullfile(ea_space,'atlases','DISTAL Minimal (Ewert 2017)','lh','STN.nii.gz'));


[I1,I2,I3] = ind2sub(STN_nii.dim, 1:numel(STN_nii.img));

q_XYZ = ea_vox2mm([I1',I2',I3'], STN_nii.mat);



STN_nii.dt = [16,0];
STN_nii.img = zeros(STN_nii.dim);

% pair-based analysis

% all recordings units  in terms of pairs)
[S_MNI_pairs,~,ix_pairs] = unique(PairLoc,'rows');
S_MNI_pairs_n = accumarray(ix_pairs,1)';
Idx_pairs = rangesearch(S_MNI_pairs, q_XYZ, res_sphere);

for type_i = 1 : numel(unit_types)
    fprintf('Preparing output for unit type %d (%d-%d) \n', unit_types(type_i),type_i,numel(unit_types))
    idxUnits = find(PairsLocation_MNI_all.S_typeFRmod == unit_types(type_i));
    [S_MNI_unittype,~,ix_unittype] = unique(PairLoc(idxUnits,:),'rows');
    S_MNI_unittype_n = accumarray(ix_unittype,1)';
    Idx_unittype = rangesearch(S_MNI_unittype, q_XYZ, res_sphere);

    TypePairdensity_STN_n = STN_nii;
    TypePairdensity_STN_n.img = zeros(STN_nii.dim);
    TypePairdensity_STN_n.img(:) = cellfun(@(x,y) sum(S_MNI_unittype_n(x))/sum(S_MNI_pairs_n(y)), Idx_unittype, Idx_pairs);
    TypePairdensity_STN_n.fname = fullfile(PATH_NIFTI,sprintf('pair_density_stn-unit-type-%s.nii',unit_types_labels{type_i}));

    if exist('ea_space') == 2
        ea_write_nii(TypePairdensity_STN_n)
    else
        spm_write_vol(TypePairdensity_STN_n,TypePairdensity_STN_n.img);
    end
    gzip(TypePairdensity_STN_n.fname);
    delete(TypePairdensity_STN_n.fname);
    fprintf('Saving nifti files \n')
end

%% cluster prob by frequency bandsd


%% cluster freq 2 again, distribution and comaprison across tasks
cluster_prob_unittypes = arrayfun(@(x) cluster_prob(ClusterS_FRMod == x, : ), unit_types,'UniformOutput',false);

nClusters_pertype = cellfun(@(x) size(x,1), cluster_prob_unittypes);
nPairs_pertype = arrayfun(@(x) sum(Unit_FRmod == x), unit_types);

nClusters_pertype_norm = nClusters_pertype./nPairs_pertype;
ClusterCentroidBand_unittypes = arrayfun(@(x) ClusterCentroidBand(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
ClusterEAtlasLabel_unittype = arrayfun(@(x) ClusterEAtlasLabel(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
uniqueClusterEAtlasLabel = unique(ClusterEAtlasLabel);

ROItoanalyze = {'G_and_S_subcentral L'
    'G_front_inf-Opercular L'
    'G_front_inf-Triangul L'
    'G_front_middle L'
    'G_pariet_inf-Supramar L'
    'G_postcentral L'
    'G_precentral L'
    'G_temp_sup-Lateral L'
    'G_temp_sup-Plan_tempo L'
    'G_temporal_middle L'};


% OLD APPROACH (this is some sort of overlap ---how many clusters together
% at the same time?)
cluster_distr_bands_unittypes = cell(1,numel(unit_types));

for type_i = 1 : numel(unit_types)
    for ti = 1 : numel(T)
        tmp = ClusterCentroidBand(cluster_prob(:,ti)>0 & ClusterS_FRMod == unit_types(type_i));
        for fi = 1:6
            cluster_distr_bands_unittypes{type_i}(fi,ti) = sum(tmp ==  fi)/nPairs_pertype(type_i);
        end


    end
end
%
% idxCluster_Task = [];
% cluster_distr_bands_task = nan(6,5);
% for evt = 1:5
%     idxCluster_Task = cat(3, T <= ClusterEvtTimes(:,evt+1) & T >= ClusterEvtTimes(:,evt));
%     cluster_distr_bands_task(:,evt) = mean(cluster_distr_bands(:,T <= ClusterEvtTimes(:,evt+1) & T >= ClusterEvtTimes(:,evt)),2,'omitnan');
% end
%idxCluster_Task = [];
% cluster_distr_task = cell(1,5);
% for evt = 1:5
%     %idxCluster_Task = cat(3, idxCluster_Task, T <= ClusterEvtTimes(:,evt+1) & T >= ClusterEvtTimes(:,evt));
%     cluster_distr_task{evt} = cluster_prob( T <= ClusterEvtTimes(:,evt+1) & T >= ClusterEvtTimes(:,evt));
% end
% cluster_distr_task = arrayfun(@(x)


cluster_distr_bands_task_unittypes = nan(4,6,5);
for type_i = 1 : numel(unit_types)
    for evt = 1:5
        cluster_distr_bands_task_unittypes(type_i,:,evt) = mean(cluster_distr_bands_unittypes{type_i}(:,T <= meanEvts(:,evt+1) & T >= meanEvts(:,evt)),2,'omitnan');
    end
end


figurer('renderer','painters')

% idxITI = T <= ClusterEvtTimes
%
% idxITI = T <= (meanCueOnset -meanSpeechOnset););
%
% idxSpeech = T <= ClusterEvtTimes(  & T  >= (meanSpeechOnset - meanSpeechOnset);
% idxCue= T <= (meanCueOffset -meanSpeechOnset)  & T  >= (meanCueOnset - meanSpeechOnset);
% idxCPreSpeech= T <= (meanSpeechOnset -meanSpeechOnset)  & T  >= (meanCueOffset - meanSpeechOnset);
% idxPostSpeech = T >= (meanSpeechOffset -meanSpeechOnset);


%
% cluster_distr_bands_task(:,1) = mean(cluster_distr_bands(:,idxITI),2,'omitnan');
% cluster_distr_bands_task(:,2) = mean(cluster_distr_bands(:,idxCue),2,'omitnan');
% cluster_distr_bands_task(:,3) = mean(cluster_distr_bands(:,idxCPreSpeech),2,'omitnan');
% cluster_distr_bands_task(:,4)= mean(cluster_distr_bands(:,idxSpeech),2,'omitnan');
% cluster_distr_bands_task(:,5) = mean(cluster_distr_bands(:,idxPostSpeech),2,'omitnan');

cluster_distr_task_untitypes = squeeze(sum(cluster_distr_bands_task_unittypes,2));

figure('Renderer','painters','Position',[300 300 400 400])
bar(1 : numel(unit_types),nClusters_pertype_norm,"FaceColor", [.8 .8 .8], 'FaceAlpha', .4);
xticks(1 : numel(unit_types))
xticklabels(unit_types_labels)
ylabel(' Clusters Occurrence [a.u] ')
box off
figname = fullfile(PATHS.saveFigures,'cluster_byunittypes');
saveas(gcf,figname,'fig')
saveas(gcf,figname,'png')
saveas(gcf,figname,'eps')


xLabeling = cfg.plot.EvtTypes(2:end-1);

figure('Renderer','painters','Position',[300 100 1200 1200])
tiledlayout(4,4)
for ii = 1 : numel(unit_types)
    nexttile(4*(ii-1) + 1,[1 2])
    hold on
    aa = area(T ,cluster_distr_bands_unittypes{ii}');
    colororder(hex2rgb(bands.color))
    %legend(aa,bands.name,'location','best')
    set(gca, 'XTick', meanEvts([2 3 4 5]))
    set(gca, 'XTickLabel', xLabeling)
    xtickangle(30)
    ylabel( {unit_types_labels{ii}, " Cluster Occurrence [a.u.] "})
    % distribution sample-by-sample frequencies
    xlim([-2.5 2])
    ylim([0 0.03])

    nexttile(4*(ii-1) + 3,[1 2])
    hold on
    plot(T ,cluster_distr_bands_unittypes{ii});
    colororder(hex2rgb(bands.color))
    %legend(aa,bands.name,'location','best')
    set(gca, 'XTick', meanEvts([2 3 4 5]))
    set(gca, 'XTickLabel', xLabeling)
    xtickangle(30)
    ylabel( " Cluster Occurrence [a.u.] ")
    % distribution sample-by-sample frequencies
    xlim([-2.5 2])
    ylim([0 0.03])

    if ii == 1
        legend(bands.name,'location','best')
    end
end

figname = fullfile(PATHS.saveFigures,'cluster_byunittypes_byfreq_bytime');
saveas(gcf,figname,'fig')
saveas(gcf,figname,'png')
saveas(gcf,figname,'eps')

%% unit type schema

ClusterSubjects_unittype = arrayfun(@(x) ClusterSubjects(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
ClusterSChannel_unittype = arrayfun(@(x) ClusterSChannel(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
SignTask_id_unittype = arrayfun(@(x) SignTask_id(:,ClusterS_FRMod == x), unit_types,'UniformOutput',false);

PairSubjects_unittype = arrayfun(@(x) PairsLocation_MNI_all.subj_id(Unit_FRmod == x), unit_types,'UniformOutput',false);
PairSChannel_unittype = arrayfun(@(x) PairsLocation_MNI_all.S_channel(Unit_FRmod == x), unit_types,'UniformOutput',false);


% figure('renderer','painters','Position',[400 100 2000 500])
% tiledlayout(1,numel(unit_types))
% for unit_i = 1 : numel(unit_types)
%     nexttile
%     imagesc(1:numel(ROItoanalyze),T, clusterROILabels_distr_bands_unittypes{type_i})
%
% end


figure('renderer','painters','Position',[400 100 2000 1500],'DefaultAxesFontSize',16)
peakband_order = cell(1,numel(unit_types));
%tiledlayout(1,numel(unit_types))
tiledlayout(2,3)
for unit_i = 1 : numel(unit_types)
    switch unit_i
        case 1
            nexttile(1,[2 1])
        case 2
            nexttile(3,[1 1])
        case 3
            nexttile(2,[2 1])
        case 4
            nexttile(6,[1 1])
    end
    hold on
    [uniqueUnits,idauniqueUnits,idxuniqueUnits] = unique([ClusterSubjects_unittype{unit_i} ClusterSChannel_unittype{unit_i}],'rows');
    nuniqueUnits_clus = size(uniqueUnits,1);
    uniqueSubjects = unique(uniqueUnits(:,1));
    nuniqueSubjects = numel(uniqueSubjects);


    nClusters_uniqueUnits = accumarray(idxuniqueUnits,1);
    modefrequencybandsUnits  = arrayfun(@(x) mode(ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == x)),unique(idxuniqueUnits));
    [modefrequencybandsUnits_sorted, idx_modefrequencybandsUnits_sorted] = sort(modefrequencybandsUnits);
    idxuniqueUnits_sorted = idxuniqueUnits(idx_modefrequencybandsUnits_sorted);
    peakband_order{unit_i} = idx_modefrequencybandsUnits_sorted;
    %uniqueUnitstrans = [1 ;find(abs(diff(idxuniqueUnits)) > 0) + 1];
    cnt = nClusters_pertype(unit_i);
    %yax = linspace(nClusters_pertype(unit_i), 0,nClusters_pertype(unit_i));
    for clust_i = 1 : nuniqueUnits_clus
        tmp_clusters = cluster_prob_unittypes{unit_i}(idxuniqueUnits == idx_modefrequencybandsUnits_sorted(clust_i),:);
        tmp_bands = ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == idx_modefrequencybandsUnits_sorted(clust_i));
        valid_clust = find(sum(tmp_clusters,2) > 0);
        plot([2.55 2.75], [cnt cnt],'color',[ 0 0 0 .4])
        %plot([2.55 2.75], [cnt cnt],'r')
        %text(2.30 ,cnt,  string(uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),1)),'FontSize',10)
        %text(2.30 ,cnt,  string(uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),2)),'FontSize',5)
        for tmp_i = valid_clust'
            plot(T(logical(tmp_clusters(tmp_i,:))), repmat(cnt,1,numel(T((logical(tmp_clusters(tmp_i,:)))))), "Color",hex2rgb(bands.color{tmp_bands(tmp_i)}),'linewidth',2.5)
            cnt = cnt - 1;
        end

    end
    axis("tight")
    %yticks(fliplr(yax(uniqueUnitstrans)));
    yticklabels(" ")
    %set(gca,'YTickLabel', flipud(ClusterUnitFR_modSubjects_unittype{1}(idauniqueUnits)),'FontSize',12)
    ylim([cnt - .5 nClusters_pertype(unit_i)+.5])

    xlim([-2.5 2])
    for evt = [2 3 4 5]'
        plot(meanEvts(evt)' *[1,1], [0 nClusters_pertype(unit_i)], '--', 'Color', 'r', 'LineWidth', 1)
    end
    set(gca, 'XTick', meanEvts([2 3 4 5]))
    xLabeling = cfg.plot.EvtTypes(2:end-1);
    %xLabeling{2} = '';
    set(gca, 'XTickLabel', xLabeling)
    xtickangle(30)
    axis tight
    ylabel([unit_types_labels{unit_i}, ' Neurons'])
    xline( 2.65)
    title(sprintf('Units %d/%d: %1.2f',nuniqueUnits_clus,nuniqueUnits_perunittype(unit_i),nClusters_pertype_norm(unit_i)))
end


% figname = fullfile(PATHS.saveFigures,'cluster_byunittypes_timeresolved');
% saveas(gcf,figname,'fig')
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'eps')



%%
SU_cluster = table();
cnt_tbl = 1;

% generate surrogates of random spectral specificity
btsp = 1000;

ClusterSingleUnitsSequence = cell(1,numel(unit_types));
check = [];
for unit_i = 1 : numel(unit_types)

    hold on
    [uniqueUnits,idauniqueUnits,idxuniqueUnits] = unique([ClusterSubjects_unittype{unit_i} ClusterSChannel_unittype{unit_i}],'rows');
    nuniqueUnits_clus = size(uniqueUnits,1);
    uniqueSubjects = unique(uniqueUnits(:,1));
    nuniqueSubjects = numel(uniqueSubjects);


    nClusters_uniqueUnits = accumarray(idxuniqueUnits,1);
    modefrequencybandsUnits  = arrayfun(@(x) mode(ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == x)),unique(idxuniqueUnits));
    [modefrequencybandsUnits_sorted, idx_modefrequencybandsUnits_sorted] = sort(modefrequencybandsUnits);
    idxuniqueUnits_sorted = idxuniqueUnits(idx_modefrequencybandsUnits_sorted);

    %uniqueUnitstrans = [1 ;find(abs(diff(idxuniqueUnits)) > 0) + 1];
    cnt = nClusters_pertype(unit_i);
    %yax = linspace(nClusters_pertype(unit_i), 0,nClusters_pertype(unit_i));
    for clust_i = idx_modefrequencybandsUnits_sorted'
        SU_cluster.SubjID(cnt_tbl) =  uniqueUnits(clust_i,1);
        SU_cluster.UnitID(cnt_tbl) =  uniqueUnits(clust_i,2);
        SU_cluster.UnitTypeFR(cnt_tbl) =  unit_types(unit_i);

        tmp_clusters = cluster_prob_unittypes{unit_i}((idxuniqueUnits == clust_i),:);
        ClusterSingleUnitsSequence{unit_i}{end+1} = tmp_clusters;
        tmp_signtask = SignTask_id_unittype{unit_i}(:,idxuniqueUnits == clust_i);
        tmp_bands = ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == clust_i);
        SU_cluster.nClusters(cnt_tbl) =  numel(tmp_bands);
        check = [check; [SU_cluster.nClusters(cnt_tbl) sum(contains(ClusterSubjects,SU_cluster.SubjID(cnt_tbl)) & strcmpi(ClusterSChannel,SU_cluster.UnitID(cnt_tbl))) ]];
        clusters_probands_byunit = arrayfun(@(x) mean(tmp_bands == x), 1:6);

        for fi = 2:6
            SU_cluster.(char(bands.name(fi)))(cnt_tbl) = clusters_probands_byunit(fi);
        end
        Entropyband = sum(-clusters_probands_byunit.*log2(clusters_probands_byunit + eps));
        if Entropyband < 0
            Entropyband = 0;
        elseif Entropyband > 1
            Entropyband = 1;
        end

        SU_cluster.BandSpecifity(cnt_tbl) = 1 - Entropyband/log2(6);
        BandSpecifity_btsp = nan(1,btsp);
        btsp_values = randi(5,btsp,numel(tmp_bands));
        btsp_values_p = nan(btsp,6);
        for btsp_i = 1 : btsp
            btsp_values_p(btsp_i,:) = arrayfun(@(x) mean(btsp_values(btsp_i,:) == x,2),1:6);
            BandSpecifity_btsp(btsp_i) = 1 - sum(-btsp_values_p(btsp_i,:).*log2(btsp_values_p(btsp_i,:) + eps))/log2(6);
        end
        SU_cluster.BandSpecifity_prct95(cnt_tbl) = prctile(BandSpecifity_btsp,95);
        SU_cluster.BandSpecifity_corr(cnt_tbl) = 1 - (Entropyband - mean(BandSpecifity_btsp,'omitnan')) /(log2(6) - mean(BandSpecifity_btsp,'omitnan')) ;
        SU_cluster.BandSpecifity_prct95_corr(cnt_tbl) = 1 - (SU_cluster.BandSpecifity_prct95(cnt_tbl)- mean(BandSpecifity_btsp,'omitnan')) /(log2(6) - mean(BandSpecifity_btsp,'omitnan')) ;

        SU_cluster.nPairs(cnt_tbl) =  sum(contains(PairSubjects_unittype{unit_i}, uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),1)) & contains(PairSChannel_unittype{unit_i},uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),2)));
        SU_cluster.ClusterDensity(cnt_tbl) =  SU_cluster.nClusters(cnt_tbl)/SU_cluster.nPairs(cnt_tbl);

        for fi = 2:6
            for evt_i = 1 : 5
                tmp = sum(tmp_signtask(:,tmp_bands == fi),2)/SU_cluster.nPairs(cnt_tbl);
                SU_cluster.([char(bands.name(fi)),'_', char(Task_labels{evt_i})])(cnt_tbl) = tmp(evt_i);
            end
        end

        % need to divide by task segment

        cnt_tbl = cnt_tbl + 1;

    end

end

SU_cluster.BandSpecifity_corr(SU_cluster.BandSpecifity_corr < 0) = 0;
SU_cluster.BandSpecifity_corr(SU_cluster.BandSpecifity_corr > 1) = 1;

% now plot the table
matrix_clusterbands = cell(1,numel(unit_types));
SU_cluster_tmp = arrayfun(@(x) SU_cluster(SU_cluster.UnitTypeFR == x,:), unit_types,'uni',false);


%%
for unit_i = 1 : numel(unit_types)
    clim_bandmaps = [0 max(max(SU_cluster_tmp{unit_i}{:,16:40}))];

    figure('Renderer','painters','Position',[300 300 1500 800])
    tiledlayout(10,7,'TileSpacing','Compact','Padding','Compact')
    nexttile(1,[9 2])
    bh = barh(SU_cluster_tmp{unit_i}{:,5:9},'stacked');
    for fi = 2:6
        bh(fi-1).FaceColor = hex2rgb(bands.color(fi));
    end
    yticks(1:height(SU_cluster_tmp{unit_i}))
    yticklabels( SU_cluster_tmp{unit_i}{:,1}   + ' ' + SU_cluster_tmp{unit_i}{:,2})
    set(gca,'YDir','Reverse')
    xlabel('Cluster proportion [a.u.]')
    box off

% 
%     nexttile(10,[9,1])
%     scatter(SU_cluster_tmp{unit_i}.BandSpecifity , 1:height(SU_cluster_tmp{unit_i}),25,'k','filled')
%     hold on
%     %hold on
%     plot(SU_cluster_tmp{unit_i}.BandSpecifity , 1:height(SU_cluster_tmp{unit_i}),'k')
%     %plot([zeros(1,height(SU_cluster_tmp{1})); SU_cluster_tmp{1}.BandSpecifity_prct95_corr'],repmat(1:height(SU_cluster_tmp{1}),2,1),'color',[0 0 0 .4],'linewidth',5 )
%     set(gca,'YDir','Reverse')
%     yticks(1:height(SU_cluster_tmp{unit_i}))
%     ylim([0 height(SU_cluster_tmp{unit_i})+1])
%     xlabel('Spectral Specificty [a.u.]')
%     yticklabels(' ')
%     xlim([0 1])

%     nexttile(11,[9,1])
%     scatter(1:height(SU_cluster_tmp{unit_i}),SU_cluster_tmp{unit_i}.BandSpecifity , 25,'k','filled')
%     hold on
%     plot(1:height(SU_cluster_tmp{unit_i}),SU_cluster_tmp{unit_i}.BandSpecifity , 'color','k','linewidth',.7)
%     set(gca,'YDir','Reverse')
%     xticks(1:height(SU_cluster_tmp{unit_i}))
%     xlim([0 height(SU_cluster_tmp{unit_i})+1])
%     ylabel('Spectral Specificty [a.u.]')
%     xticklabels(' ')
%     ylim([0 1])
%         box off
% 


%      
%     scatter(SU_cluster_tmp{unit_i}.ClusterDensity , 1:height(SU_cluster_tmp{unit_i}),25,'k','filled')
% ax2.XAxisLocation = 'top';
%     plot( 1:height(SU_cluster_tmp{unit_i}),SU_cluster_tmp{unit_i}.ClusterDensity ,'k')
%     set(gca,'YDir','Reverse')
%     xlabel('Cluster Density [a.u.]')
%     yticklabels(' ')
%     xlim([0 1])
%     box off
    %=camroll(90)

%     nexttile(1,[1 2])
%     bh = barh([mean(SU_cluster_tmp{unit_i}{:,5:9}); nan(1,5)],'stacked');
%     for fi = 2:6
%         bh(fi-1).FaceColor = hex2rgb(bands.color(fi));
%         bh(fi-1).BarWidth = 115;
%     end
%     ylim([0.5 1.5])
%     box off
%     axis off

    for fi = 2:6
        
        nexttile(3 + (fi -2),[9 1])
        %figure
        imagesc(1:5,1:height(SU_cluster_tmp{unit_i}),SU_cluster_tmp{unit_i}{:,{char(strcat(bands.name(fi),'_ITI')),char(strcat(bands.name(fi),'_Cue')), ...
            char(strcat(bands.name(fi),'_Pre-Speech')), char(strcat(bands.name(fi),'_Speech')), char(strcat(bands.name(fi),'_Post-Speech'))}});
        ylim([0 height(SU_cluster_tmp{unit_i})+1])

        xticks(1:5)
        set(gca, 'XTickLabel', Task_labels)
        xtickangle(30)
        yticks(1:height(SU_cluster_tmp{unit_i}))
        yticklabels(' ')
        colormap(flipud(gray))
        ylim([0 height(SU_cluster_tmp{unit_i})+1])
        set(gca,'TickDir','in')
        clim(clim_bandmaps)
        title(bands.symbol(fi),'FontSize',14)
    end
    cb = colorbar;
    ylabel(cb,'Cluster Density [a.u.]','FontSize',10)
    

    figname = fullfile(PATHS.saveFigures,['cluster_byunittypes_unitresolved_',char(unit_types_labels(unit_i))]);
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'png')
    saveas(gcf,figname,'eps')

end

%%

clusterrecap_mean = nan(numel(unit_types),5,5);
clusterrecap_sem = nan(numel(unit_types),5,5);

for unit_i = 1 : numel(unit_types)
    for fi = 2:6
        clusterrecap_mean(unit_i,fi-1,:) = mean(SU_cluster_tmp{unit_i}{:,{char(strcat(bands.name(fi),'_ITI')),char(strcat(bands.name(fi),'_Cue')), ...
            char(strcat(bands.name(fi),'_Pre-Speech')), char(strcat(bands.name(fi),'_Speech')), char(strcat(bands.name(fi),'_Post-Speech'))}});
        clusterrecap_sem(unit_i,fi-1,:) = std(SU_cluster_tmp{unit_i}{:,{char(strcat(bands.name(fi),'_ITI')),char(strcat(bands.name(fi),'_Cue')), ...
            char(strcat(bands.name(fi),'_Pre-Speech')), char(strcat(bands.name(fi),'_Speech')), char(strcat(bands.name(fi),'_Post-Speech'))}})/sqrt(numel(SU_cluster_tmp{unit_i}));
    end
end

figure('renderer','painters','position',[300 300 1400 400])
for unit_i = 1:4
    nexttile
    hold on
    for fi = 2:6
        plot(squeeze(clusterrecap_mean(unit_i,fi-1,:))','color',hex2rgb(bands.color(fi)))
    end
end

% 
% if size(y,1)==2 %plot shaded area
%     px=[x,fliplr(x)]; % make closed patch
%     py=[y(1,:), fliplr(y(2,:))];
%     patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
% end;
%  
% if size(y,1)==3 % also draw mean
%     px=[x,fliplr(x)];
%     py=[y(1,:), fliplr(y(3,:))];
%     patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
%     plot(x,y(2,:),fstr);
% end;
%%
figure('renderer','painters')
plim = 0.4;
%
scatter(SU_cluster.beta_ITI,SU_cluster.("beta_Post-Speech"),35,'k','filled')
xlim([0 plim-0.1])
ylim([0 plim])
hold on
plot([0 plim], [0 plim],'--k')

figure('renderer','painters')
plim = 0.4;
histogram(SU_cluster.beta_ITI - SU_cluster.("beta_Post-Speech"),25)
xlim([-1 1])
