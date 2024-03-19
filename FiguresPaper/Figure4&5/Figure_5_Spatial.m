  
clc
clear all
close all

%% Step 0: settings of the script
Figure_4_settings;

%% Step 1. gather all SPC (pairs & clusters) information
DB_all = get_SPCdata(SUBJECTS,'all', PATHS);

%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_all.Pairs);
% UnitOfInt = [-1 0 1 2];
% UnitLabel = ['D','N','I','M'];
% 3.1a Show PPCz on Speech
clust_res = .005;
T = [-2.5 2];

%% Step 3 Get Nifti files (nodz and nifti) 
targets = struct();
targets.Beta_Horn2017 = [-14.4, -13.2, -4.8];
targets.Alpha_Horn2017 = [-13.5,-10.3,-2.9];
targets.DBSOpt_Caire2013 = [-12.6 -13.4 -5.9];
targets.ThetaAlphaCoh_Wijk2022 = [-12.5 -16.35 -7.66];
targets.BetaCoh_Wijk2022 = [-11.85 -14.18 -6];

cfg_sp = struct();
cfg_sp.bands = bands;
cfg_sp.map = fullfile(PATH_OUTPUT,'Figures','default','group_freq','AllFrequencies','STN','STN_coords','STN_SPCdensity_roi-mm_time-none.csv');
cfg_sp.ref = fullfile(PATH_DISTAL,'atlas_index.mat');
cfg_sp.targets = targets;
cfg_sp.nperms = 500;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;
cfg_sp.EvtTimes = EvtTimes.Speech;
cfg_sp.T = T(1) : clust_res : T(2);
cfg_sp.cfg = cfg;

fh = create_SPCSTN_CSFS(cfg_sp, DB_all);


%%
for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,['SPC_Clusters_ECoG_roi-',ROI_AtlasLabels{fh_i}]);
    saveFigures(fh{fh_i}, figname)
end
clear fh


%% check if SPC maps overlap with DBS target
%cfg_sp.ref = fullfile(PATH_DISTAL,'atlas_index.mat');
%cfg_sp.atlas = 'E_atlas_label_Destrieux';

cfg = [];
cfg.VarFields = {'ClusterCentroidBand','ClusterS_MNI'};
[ClusterCentroidBand,ClusterS_MNI] = get_SPCClustersProperties(cfg,DB_all.Clusters);
nperms = 500;
PairsLoc = DB_all.Pairs.Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}};

STN_spots = table('Size',[6 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',{'X','Y','Z','Color','Size'});

for fi = 1:5
    coords = ClusterS_MNI(ClusterCentroidBand == fi+1,:);
    STN_spots{fi,{'X','Y','Z'}} = mean(coords,'omitnan');
    STN_spots.Size(fi) = 1;
    STN_spots.Color(fi) = fi;
end
STN_spots{fi + 1,{'X','Y','Z'}} = targets.DBSOpt_Caire2013;
STN_spots.Size(fi+1) = 2;
STN_spots.Color(fi+1) = fi+1;