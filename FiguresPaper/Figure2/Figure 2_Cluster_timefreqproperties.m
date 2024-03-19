clc
clear all
close all

%% Step 0: settings of the script
Figure_2_settings;

%% Step 1. gather all SPC (pairs & clusters) information
DB_all = get_SPCdata(SUBJECTS,'all', PATHS);

%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_all.Pairs);
clust_res = .005;
T = -2.5 : clust_res : 2;
%% Step 3. Get Cluster Information that I need


%[PPC,PPCz,clust_prob] = get_SPCClustersProperties(cfg_cl, DB_all.Clusters);



%% Step 4. Time & frequency properties

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.VarFields = {'ClusterCentroid','ClusterCentroidBand','ClusterOnOff'};
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;

% plot %nclusters x time
fh = plot_ClusterTime(cfg_cl,DB_all.Clusters);

figname = fullfile(PATHS.saveFigures,'SPC_Clusters_timefreq');
saveFigures(fh, figname)

%% Step 5: Plot Duration and Frequency
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'Duration_1DHist','Frequency_1DHist'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);

for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,['SPC_',cfg_cl.Properties{fh_i}]);
    saveFigures(fh{fh_i}, figname)
end


%% Step 6: Cluster frequency time distribution

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DArea'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);
figname = fullfile(PATHS.saveFigures,['SPC_','TaskModulation_1DArea']);
saveFigures(fh, figname)

%% Step 7: Cluster robustness
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'RobustnessCheck'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);
figname = fullfile(PATHS.saveFigures,['SPC_','RobustnessCheck']);
saveFigures(fh{1}, figname)

%% Step 8: Cluster aggregation
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DAggregationLine'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);

figname = fullfile(PATHS.saveFigures,['SPC_','TaskModulation_1DAggregationLine']);
saveFigures(fh{1}, figname)


%% Step 9: Cluster phase population agreement
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DPhaseAggregation'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);

figname = fullfile(PATHS.saveFigures,['SPC_','TaskModulation_1DPhaseAggregation']);
saveFigures(fh{1}, figname)

%% Step 10: Clusters lag delay
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.freq_range = [4 12];
cfg_cl.Properties = {'TaskModulation_LagDelay','TaskModulation_Phasepolar','TaskModulation_Lagdelay_overtime'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);

for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,['SPC_low',cfg_cl.Properties{fh_i}]);
    saveFigures(fh{fh_i}, figname)
end



cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.freq_range = [4 30];
cfg_cl.Properties = {'TaskModulation_LagDelay','TaskModulation_Phasepolar','TaskModulation_Lagdelay_overtime'};
fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters);

for fh_i = 1 : numel(fh)
    figname = fullfile(PATHS.saveFigures,['SPC_beta',cfg_cl.Properties{fh_i}]);
    saveFigures(fh{fh_i}, figname)
end
%% Step 11: Cluster lag delay per time