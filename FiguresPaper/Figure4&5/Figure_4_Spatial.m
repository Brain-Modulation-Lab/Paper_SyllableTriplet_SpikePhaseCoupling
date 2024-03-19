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

cfg_sp = struct();
cfg_sp.time = 'all'; 
cfg_sp.cfg = cfg;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;
cfg_sp.clust_res = clust_res;
cfg_sp.nperms = 500;
cfg_sp.bands = bands;
cfg_sp.EvtTimes = EvtTimes.Speech;
cfg_sp.T = T(1) : clust_res : T(2);
[fh, ROI_AtlasLabels] = create_SPCROI_time(cfg_sp, DB_all);

for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,['SPC_Clusters_ECoG_roi-',ROI_AtlasLabels{fh_i}]);
    saveFigures(fh{fh_i}, figname)
end
clear fh

%% try plot density w.r.t. distances

cfg_sp = struct();
cfg_sp.bands = bands;
cfg_sp.map = fullfile(PATH_OUTPUT,'Figures','default','group_freq','res-1mm','nii','clus_density_res-1mm_band-group_timeframe-avg_reduced.txt');
cfg_sp.ref = fullfile(PATH_UTILITIES,'resources','cortex_MNI.mat');
cfg_sp.nperms = 500;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;
fh = create_SPCCortex_CSFS(cfg_sp, DB_all);

fignames = {'SPC_ROI_2Dviz','SPC_ROI_density','SPC_ROI_duration','SPC_ROI_freq','SPC_ROI_focality','SPC_ROI_UnitFT_density_simple','SPC_ROI_UnitFT_density_refined','SPC_ROI_densitypermutation_bands'}; 
for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,fignames{fh_i});
    saveFigures(fh{fh_i}, figname)
end
clear fh

%% try A-P and D-V directions
cfg_sp = struct();
cfg_sp.bands = bands;
cfg_sp.map = fullfile(PATH_OUTPUT,'Figures','default','group_freq','res-1mm','nii','clus_density_res-1mm_band-group_timeframe-avg_reduced.txt');
cfg_sp.ref = fullfile(PATH_UTILITIES,'resources','cortex_MNI.mat');
cfg_sp.nperms = 500;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;

fh = create_SPCDensity_distr(cfg_sp,DB_all);
figname = fullfile(PATHS.saveFigures,'SPCDensity_ECOG_planes');

saveFigures(fh, figname)


%% check comparison with GODIVA
cfg_diva = struct();
cfg_diva.map = fullfile(PATH_OUTPUT,'Figures','default','group_freq','AllFrequencies','Cortex','Cortex_2mm_nii','clus_density_res-1mm_band-group_timeframe-avg_reduced.txt');
cfg_diva.ref_coords = 'DIVA_locationmaps_coverage.xlsx';
cfg_diva.bands = bands;
cfg_diva.nperms = 500;

fh = compare_DIVAlocation(cfg_diva, DB_all);
figname = fullfile(PATHS.saveFigures,'SPC_DIVA_comparison');
saveFigures(fh, figname)


