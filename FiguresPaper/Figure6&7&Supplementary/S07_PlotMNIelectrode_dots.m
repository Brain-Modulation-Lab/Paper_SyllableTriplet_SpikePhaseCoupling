clc
clear all
close all


if ismac
    PATH_SERVER = '/Volumes/Nexus/DBS';
    PATH_PROJECT = '/Volumes/Nexus4/Users/MV1019/PhaseLocking';
    PATH_RESOURCES = '/Volumes/Nexus1/Resources/MNI_Cortex_plotting';
elseif ispc
    PATH_SERVER = 'Z:\DBS';
    PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
    PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
end


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

cfg = set_configs("default");
TYPE_DB = 'default';
fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

PATHS = struct();
PATHS.DataDir      = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures_new',TYPE_DB);



bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
%% Step 1. gather all clusters in a very big structure
toDo = [1 : n_SUBJECTS];% 1 : n_SUBJECTS;
[nSign_pairs_all, nPairs_all]  = deal(nan(1,n_SUBJECTS)) ;
Clusters_all = [];
PPC_mat_all = [];
PPCz_mat_all = [];
PPCmedianperm_mat_all = [];
ES_mat_all = [];
Phase_mat_all = [];
PairsLocation_MNI_all = [];
PLVTimeEvts_all = [];
toCorr  = 0;
for subj_i = toDo %
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    
    PATH_clusters = fullfile(PATHS.saveMatFiles,SUBJECT,'ClustersPLV.mat');
    
    if isfile(PATH_clusters)
        fprintf(' Stacking clusters results in %s \n', PATH_clusters)
        load(PATH_clusters,'Clusters', 'nSign_pairs','n_pairs','PairsLocation_MNI', 'PPC_mat','PPCmedianperm_mat','PPCz_mat','ES_mat','PLVTimeEvts','Phase_mat');
        % put info about subjects id
        [Clusters.subj_id] = deal(SUBJECT);
        if subj_i  > 1 && subj_i ~= 6
            Clusters = arrayfun(@(x) setfield(x,'pair_i_all',x.pair_i +  toCorr),Clusters);
        elseif subj_i == 6
            Clusters = arrayfun(@(x) setfield(x,'pair_i_all',x.pair_i +  toCorr),Clusters);
        else
            [Clusters.pair_i_all] = deal(Clusters.pair_i);
        end
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
        toCorr =toCorr + n_pairs;
            
        fprintf(' Completed Stacking clusters results in %s \n', PATH_clusters)
        % clear variables for consistency
        clear Clusters nSign_pairs n_pairs PairsLocation_MNI PPC_mat PPCz_mat ES_mat PPCmedianperm_mat PLVTimeEvts Phase_mat
    else
        warning('Analysis still running: Clusters  is not available yet for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    end
end

nClusters = numel(Clusters_all);
nonan = ~isnan(nSign_pairs_all); % to manage the lack of data for patient 5

%%

E_MNI = [PairsLocation_MNI_all.E_MNI_X PairsLocation_MNI_all.E_MNI_Y PairsLocation_MNI_all.E_MNI_Z];
E_MNI_unique = unique(E_MNI,'rows');


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
ClustersPairs = [];


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
    ClustersPairs = [ClustersPairs; repmat(Clusters_all(clus_i).pair_i_all,ClusternSign(clus_i),1)];
    
end
ClusterSubjects = cellstr(ClusterSubjects);
nClusters_all = size(ClusterOnOff,1);

ClusterCentroidBand = nan(1,nClusters_all);
for clus_i = 1 : nClusters_all
    ClusterCentroidBand(clus_i) = find(ClusterCentroid(clus_i,2) < bands.fends,1);
end