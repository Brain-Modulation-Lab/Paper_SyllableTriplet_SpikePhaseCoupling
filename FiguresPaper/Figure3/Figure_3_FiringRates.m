clc
clear all
close all

%% Step 0: settings of the script
Figure_3_settings;

%% Step 1. gather all SPC (pairs & clusters) information
DB_all = get_SPCdata(SUBJECTS,'all', PATHS);

%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_all.Pairs);
UnitOfInt = [-1 0 1 2];
UnitLabel = ['D','N','I','M'];
% 3.1a Show PPCz on Speech
clust_res = .005;
T = [-2.5 2];
%% Step 3. Get information of firing rate modulations

%% Step 4. Calculate maps for each units type

% Increasing Type (Unit_Info == 1);
% None Type (Unit_Info == 0);
% Decreasing Type (Unit_Info == -1);
% Mixed Type (Unit_Info == 2);

R
nPerms = 0;

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-.02 1.8];
cfg_fh.CMap = linspecer;
cfg_fh.NClustplot = 'all';
cfg_fh.Baseline = [-Inf EvtTimes.Speech(2)];



for unit_type = 1 : numel(UnitOfInt)
    PairType = [UnitLabel(unit_type),'OnlySign'];
    SPCmap = get_SPCmap(DB_all.Pairs,'PPCz', T, ...
        Tres,cfg,'Speech',PairType,'default',nPerms);
    fh = plotter_SPCmap(SPCmap, cfg_fh);
    if nPerms > 0 && numel(fh) == 2
        figname = fullfile(PATHS.saveFigures,['PPCzmap_pairs-unit-', PairType, '_evt-speech']);
        saveFigures(fh{1}, figname)
        figname = fullfile(PATHS.saveFigures,['PPCtstat_pairs-unit-', PairType, '_evt-speech']);
        saveFigures(fh{2}, figname)
    else
        figname = fullfile(PATHS.saveFigures,['PPCzmap_pairs-unit-', PairType, '_evt-speech']);
        try
            saveFigures(fh, figname)
        catch
            saveFigures(fh{1}, figname)
        end
    end
end

%% Step 5: Calculate Cluster time x frequency

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.VarFields = {'ClusterCentroid','ClusterCentroidBand','ClusterOnOff','ClusterS_FRMod'};
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.unit_type = UnitLabel;
cfg_cl.order = 'freq-wise';


%Unit_Info = get_FRTypes(DB_all.Pairs);
% plot %nclusters x time
fh = plot_ClusterTime(cfg_cl,DB_all.Clusters);
for unit_type = 1 : numel(UnitOfInt)
        figname = fullfile(PATHS.saveFigures,['SPC_Clusters-unit-', UnitLabel(unit_type), '_timefreq']);
        saveFigures(fh{unit_type}, figname)
end


%% Step 6: Plot Duration and Frequency
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'Duration_1DHist','Frequency_1DHist'};
%cfg_cl.VarFields = {'ClusterS_FRMod'};

%ClusterS_FRMod = get_SPCClustersProperties(cfg_cl, DB_all.Clusters);
for unit_type = 1 : numel(UnitOfInt)
    fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters([DB_all.Clusters.S_typeFRmod] == UnitOfInt(unit_type)));
    for fh_i = 1 : numel(fh)
        figname = fullfile(PATHS.saveFigures,['SPC_Clusters-unit-',UnitLabel(unit_type),'_',cfg_cl.Properties{fh_i}]);
        saveFigures(fh{fh_i}, figname)
    end
    clear fh
end

%% Step 7: Cluster frequency time distribution

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T(1) : Tres : T(2);
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DAggregationLine'};

for unit_type = 1 : numel(UnitOfInt)
    fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters([DB_all.Clusters.S_typeFRmod] == UnitOfInt(unit_type)));
    for fh_i = 1 : numel(fh)
        figname = fullfile(PATHS.saveFigures,['SPC_Clusters-unit-',UnitLabel(unit_type),'_',cfg_cl.Properties{fh_i}]);
        try
            saveFigures(fh, figname)
        catch
            saveFigures(fh{1}, figname)
        end
    end
    clear fh
end

%% Step 8: plot lag (this is a try)
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T(1) : Tres : T(2);
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_LagDelay','TaskModulation_Phasepolar','TaskModulation_Lagdelay_overtime'};
cfg_cl.freq_range = [12 30];

for unit_type = 1 : numel(UnitOfInt)
    fh = plot_ClusterProperties(cfg_cl, DB_all.Clusters([DB_all.Clusters.S_typeFRmod] == UnitOfInt(unit_type)));
    for fh_i = 1 : numel(fh)
        figname = fullfile(PATHS.saveFigures,['SPC_',cfg_cl.Properties{fh_i}]);
        try
            saveFigures(fh, figname)
        catch
            saveFigures(fh{1}, figname)
        end
    end
    clear fh
end

%% plot occurrence by cell types
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T(1) : Tres : T(2);
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;

PairType = arrayfun(@(x) DB_all.Pairs.Location([DB_all.Pairs.Location.S_typeFRmod] == x,:), UnitOfInt,'uni',false);
ClustType = arrayfun(@(x) DB_all.Clusters([DB_all.Clusters.S_typeFRmod] == x), UnitOfInt,'uni',false);

ClustDensityType = cellfun(@(x,y) sum([x.nSign])/height(y)*100,ClustType,PairType);
PairSignType = cellfun(@(x) mean(x.nSign > 0),PairType);

cfg_cl.VarFields = {'ClusterS_FRMod','ClusterCentroidBand'};
[ClusterS_FRMod,ClusterCentroidBand] = get_SPCClustersProperties(cfg_cl,DB_all.Clusters);

ClustDensityTypeByBand = nan(4,6);
for ui = 1:4
    for fi = 2:6
        ClustDensityTypeByBand(ui,fi) = sum(ClusterS_FRMod' == UnitOfInt(ui) & ClusterCentroidBand == fi)/height(PairType{ui})*100;
    end
end

fh = figure('position',[200 200 400 600]);
tiledlayout(4,1)
nexttile(1,[1 1])
    imagesc(ClustDensityType )
    box off
    xticks(1:4);
    yticks(1);
    yticklabels({'All'});
    xticklabels({'D','N','I','M'}  )
    colormap(linspecer);
    cb = colorbar;
    cb.Label.String = 't-SPC density';
    nexttile(2,[3,1])
    imagesc(ClustDensityTypeByBand(:,2:6)' )
    box off
    xticks(1:4);
    yticks(1:5);
    yticklabels(bands.symbol(2:end));
    xticklabels({'D','N','I','M'}  )
    colormap(linspecer);
    cb = colorbar;
    cb.Label.String = 't-SPC density';


% for ui = 1:4
%     nexttile
%     bar(1:6,[ClustDensityType(ui) ClustDensityTypeByBand(ui,2:6)])
%     box off
%     xticks(1:6);
%     xticklabels({'All','theta','alpha','beta','gammaL','gammaH'})
%     ylabel('t-SPC density')
%     ylim([0 25])
% end
figname = fullfile(PATHS.saveFigures,['SPC_Density_UnitType']);
saveFigures(fh, figname)


%ClustDensityTypeByType = cellfun(@(x,y) sum([x.nSign])/height(y)*100,ClustType,PairType);

%% plot baseline vs iti cells

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T(1) : Tres : T(2);
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
%cfg_cl.Properties = {'TaskModulation_LagDelay','TaskModulation_Phasepolar','TaskModulation_Lagdelay_overtime'};
%cfg_cl.freq_range = [12 30];

fh = plot_SingleUnits(cfg_cl,DB_all);

figname = fullfile(PATHS.saveFigures,['SPC_Clusters-unit-',UnitLabel(unit_type),'_',cfg_cl.Properties{fh_i}]);
try
    saveFigures(fh, figname)
catch
    saveFigures(fh{1}, figname)
end
figname = fullfile(PATHS.saveFigures,['SPC_Clusters_biplot-unit-',UnitLabel(unit_type),'_',cfg_cl.Properties{fh_i}]);
saveFigures(fh{2}, figname)

%% extract betaprofile for each cell

% get units pairs
[Units, idxUnitPairs] = getUniquePairs(DB_all.Pairs.Location);
nUnits = size(Units,1);
nUniquePairs = accumarray(idxUnitPairs,1);

% get unit clusters
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = T(1) : Tres : T(2);
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.VarFields = {'Sign_Taskid','ClusterCentroidBand','ClusterSChannel','ClusterSubjects','ClusterS_FRMod','ClusterOnOff','ClusterCentroid','clust_prob'};
[Sign_Task,ClusterCentroidBand,Cluster_SUnit, Cluster_Subject, ClusterS_FRMod,ClusterOnOff,ClusterCentroid,cluster_prob] = get_SPCClustersProperties(cfg_cl,DB_all.Clusters);
%Units_FRmod = arrayfun(@(x) , unique(idxUnitPairs));

%% check relationship with firing rate

load('FRDATA_group.mat');
SPCDensityUnits = readtable(fullfile(PATHS.saveMatFiles,"SPCDensity_units.txt"));
UnitsFRProp = strings(numel(DATA),1);
UnitID = [DATA.unit_id]';
FRbase = nan(1,numel(DATA));
FRmean = nan(1,numel(DATA));

FRzscore = nan(1,numel(DATA));

for ui = 1 : numel(DATA)
    UnitsFRProp(ui) = strcat("Unit",strrep(num2str(UnitID(ui)),' ',''));

    
      FRbase(ui) = DATA(ui).SpeechOnset.basemeanIFR;
      FRmean(ui) = mean(DATA(ui).SpeechOnset.meanifr,'omitnan');
      time = linspace(DATA(ui).SpeechOnset.respInterval(1),DATA(ui).SpeechOnset.respInterval(2),size(DATA(ui).SpeechOnset.DD,2));

      FRzscoretime = (DATA(ui).SpeechOnset.meanifr -  FRbase(ui))/std(mean(DATA(ui).SpeechOnset.IFRbase,1),[],2);
      FRzscore(ui) = max(FRzscoretime(time >=0 & time <= 1.4));
%     FRzscoretime = (DATA(ui).SpeechOnset.meanifr - mean(DATA(ui).SpeechOnset.IFRbase(:),'omitnan'))/std(DATA(ui).SpeechOnset.IFRbase(:),[],'omitnan');
%     [~,idxchange_max] = max(abs(FRchangetime));
%     [~,idxzscore_max] = max(abs(FRzscoretime));
%     FRchange(ui) = FRchangetime(idxchange_max);
%     FRchange_abs(ui) = abs(FRchangetime(idxchange_max));
%     FRzscore(ui) = FRzscoretime(idxzscore_max);
end
UnitsFRProp = [UnitsFRProp {DATA.SubjectID}'];

for xi = 1 : height(SPCDensityUnits)
    [~,idx_] =ismember(UnitsFRProp,[SPCDensityUnits.Unit(xi) SPCDensityUnits.Subject(xi)],'rows');
    try
        SPCDensityUnits.FRbase(xi) = FRbase(find(idx_));
        SPCDensityUnits.FRmean(xi) = FRmean(find(idx_));
%         SPCDensityUnits.FRchange(xi) = FRchange(find(idx_));
%         SPCDensityUnits.FRchange_abs(xi) = FRchange_abs(find(idx_));
        SPCDensityUnits.FRzscore(xi) = FRzscore(find(idx_));
    catch
        SPCDensityUnits.FRbase(xi) = nan;
        SPCDensityUnits.FRmean(xi) = nan;
%         SPCDensityUnits.FRchange(xi) = nan;
%         SPCDensityUnits.FRchange_abs(xi) = nan;
        SPCDensityUnits.FRzscore(xi) = nan;
    end
end

FR_basetype = arrayfun(@(x) FRbase(SPCDensityUnits.FRmod == x),[-1 0 1 2],'uni',false);
FR_meantype = arrayfun(@(x) FRmean(SPCDensityUnits.FRmod == x),[-1 0 1 2],'uni',false);

cellfun(@(x) mean(x),FR_basetype);
cellfun(@(x) std(x),FR_basetype);
cellfun(@(x) mean(x),FR_meantype);
cellfun(@(x) std(x),FR_meantype);

BetaDesynch = SPCDensityUnits.beta_speech - SPCDensityUnits.beta_ITI;
ThetaAlphaSynch = SPCDensityUnits.thetaalpha_speech - SPCDensityUnits.thetaalpha_ITI;
BetaRebound = SPCDensityUnits.beta_postspeech - SPCDensityUnits.beta_ITI;

% compare all units types together doing plots

cluster_probBeta_units = zeros(nUnits, numel(cfg_cl.T));
iti_speech = zeros(nUnits,4);

for unit_i = 1 : nUnits
    idxU = all(contains([Cluster_SUnit string(Cluster_Subject)],Units(unit_i,:)),2);
    betaidxU = ClusterCentroidBand' == 4 & idxU;
    if sum(betaidxU) > 0
        cluster_probBeta_units(unit_i,:) = sum(cluster_prob(betaidxU,:),1)/nUniquePairs(unit_i)*100;
    end
end
figure('position',[200 200 800 800])
tiledlayout(2,2)
nexttile(1,[1, 2])
plot(cfg_cl.T,cluster_probBeta_units' ,'k')
xlabel(' Time [s] ')
ylabel('t-SPC density ')
xlim([cfg_cl.T(1) cfg_cl.T(end)])
box off

full_betadesynch =  SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0;
hold on
plot(cfg_cl.T,cluster_probBeta_units(full_betadesynch,:)' ,'r')



% beta ITI vs speech betarange
offset = 0.01;

nexttile
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset), log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
% hold on
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1)+ offset),75,'^k','filled')
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,log10( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2)+ offset),75,'ok','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset), (SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
hold on
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , (SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , (SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1)+ offset),75,'^k','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2)+ offset),75,'ok','filled')
% highlight desynch cell
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'vr','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'xr')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'^r','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'or','filled')

xlabel('beta SPC ITI')
ylabel('beta SPC speech')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.005 100])
ylim([0.005 100])
plot([0.005 100],[0.005 100],'--k')
% xlim([-17 2])
% ylim([-17 2])
% plot([-17 2],[-17 2],'--k')

% bneta ITI vs beta rebound
offset = 0.01;
nexttile
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset), log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
% hold on
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1)+ offset),75,'^k','filled')
% scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,log10( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2)+ offset),75,'ok','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset), (SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
hold on
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , (SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , (SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 1)+ offset),75,'^k','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 2)+ offset),75,'ok','filled')
% highlight desynch cell
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'vr','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'xr')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'^r','filled')
scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'or','filled')

xlabel('beta SPC ITI')
ylabel('beta SPC postspeech')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.005 100])
ylim([0.005 100])
plot([0.005 100],[0.005 100],'--k')


% convert desynch vs rebound
% nexttile
% % bneta ITI vs beta rebound
% offset = 0.01;
% % scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset), log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
% % hold on
% % scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
% % scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1)+ offset),75,'^k','filled')
% % scatter(log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,log10( SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2)+ offset),75,'ok','filled')
% scatter(log10((SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1)+ offset))/(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset)), log10(SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == -1) )- log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) + offset),75,'vk','filled')
% hold on
% scatter(log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0)+ offset ) - log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0)+ offset) , log10(SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 0))- log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0) + offset),75,'xk')
% scatter(log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1)+ offset ) - log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1)+ offset) , log10(SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 1))- log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1) + offset),75,'^k','filled')
% scatter(log10(SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2)+ offset ) - log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2)+ offset) ,log10( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 2))- log10(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2) + offset),75,'ok','filled')
% % % highlight desynch cell
% % scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == -1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'vr','filled')
% % scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 0 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'xr')
% % scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 1 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'^r','filled')
% % scatter((SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset) ,( SPCDensityUnits.beta_postspeech(SPCDensityUnits.FRmod == 2 & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0)+ offset),75,'or','filled')
% 
% xlabel('beta SPC synch speech')
% ylabel('beta SPC synch postspeech')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlim([0.005 100])
% ylim([0.005 100])
% plot([0.005 100],[0.005 100],'--k')





figure('position',[200 200 600 600])
scatter(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == -1) , SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == -1) ,75,'vk','filled')
hold on
scatter(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 0), SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 0) ,75,'xk')
scatter(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 1) , SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 1),75,'^k','filled')
scatter(SPCDensityUnits.beta_ITI(SPCDensityUnits.FRmod == 2) ,SPCDensityUnits.beta_speech(SPCDensityUnits.FRmod == 2),75,'ok','filled')

xlabel('beta SPC ITI')
ylabel('beta SPC speech')
xlim([-17 2])
ylim([-17 2])
plot([-17 2],[-17 2],'--k')



% figure('position',[200 200 600 600])
% scatter(log10(SPCDensityUnits.beta_ITI + eps), log10(SPCDensityUnits.beta_speech + eps))
% set(gca,'xscale','log')
%  set(gca,'yscale','log')
%  xlim([10^-17 1])

%%
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRbase)),SPCDensityUnits.FRbase(~isnan(SPCDensityUnits.FRbase)),'type','Spearman')
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRbase)),SPCDensityUnits.FRmean(~isnan(SPCDensityUnits.FRbase)),'type','Spearman')
% 
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRmean)),SPCDensityUnits.FRmean(~isnan(SPCDensityUnits.FRmean)),'type','Spearman')
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRchange)),SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRchange)),'type','Spearman')
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRchange_abs)),SPCDensityUnits.FRchange_abs(~isnan(SPCDensityUnits.FRchange_abs)),'type','Spearman')
% [r,p] = corr(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRzscore)),SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRzscore)),'type','Spearman')

[r,~,p_perm,fh] = corr_permute(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRmean)),SPCDensityUnits.FRbase(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);

axfh = fh.CurrentAxes;
figure('Position',[200 200 800 300])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRmean)),SPCDensityUnits.FRmean(~isnan(SPCDensityUnits.FRmean)),'k','filled')
xlabel('t-SPC density')
ylabel('IFR [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

% compare frequency
[r,~,p_perm,fh] = corr_permute(SPCDensityUnits.meanfreq(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),SPCDensityUnits.FRbase(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),"Spearman",500,0);
%[r,~,p_perm,fh] = corr_permute(SPCDensityUnits.medianfreq(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.medianfreq)),SPCDensityUnits.FRmean(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.medianfreq)),"Spearman",500,0)

axfh = fh.CurrentAxes;

figure('Position',[200 200 800 300])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.meanfreq(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),SPCDensityUnits.FRmean(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),'k','filled')
xlabel('t-SPC frequency [Hz]')
ylabel('IFR [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;


[r,~,p_perm,fh] = corr_permute(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRmean)),SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);

axfh = fh.CurrentAxes;

figure('Position',[200 200 800 300])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.Density(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & ~isnan(SPCDensityUnits.meanfreq)),'k','filled')
xlabel('t-SPC frequency [Hz]')
ylabel('IFR [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

% check within neurons
% Tvec = T(1) : clust_res : T(2);
% cfg.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
% cfg.VarFields = {'clust_prob','Sign_Taskid','ClusterCentroidBand'};
% cfg.T = Tvec;
% [clust_prob,Sign_Taskid,ClusterCentroidBand] = get_SPCClustersProperties(cfg,DB_all.Clusters);
% 
% [Units, idxUnits] = getUniquePairs(DB_all.Clusters);
% cluster_probunits = nan(height(SPCDensityUnits),numel(Tvec));
% cluster_probunits_window = zeros(height(SPCDensityUnits),5);
% cluster_probunits_window_nogamma = zeros(height(SPCDensityUnits),5);
% cluster_probunits_window_thetaalpha = zeros(height(SPCDensityUnits),5);
% cluster_probunits_window_beta = zeros(height(SPCDensityUnits),5);
% 
% for xi = 1 : height(SPCDensityUnits)
%     [~,idx_] = ismember(Units,[SPCDensityUnits.Subject(xi) SPCDensityUnits.Unit(xi)],'rows');
%     idx_ = find(idx_);
%     if ~isempty(idx_)
%         cluster_probunits(xi,:) = 100*mean(clust_prob(idxUnits == idx_,:),1);
%         for evt_i = 1:5
%             cluster_probunits_window(xi,evt_i) = sum(Sign_Taskid(evt_i,idxUnits == idx_),2)*100/SPCDensityUnits.nPairs(xi);
%             if sum(idxUnits == idx_ & ClusterCentroidBand' < 5) > 0
%             cluster_probunits_window_nogamma(xi,evt_i) = sum(Sign_Taskid(evt_i,idxUnits == idx_ & ClusterCentroidBand' < 5),2)*100/SPCDensityUnits.nPairs(xi);
%             end
%             if sum(idxUnits == idx_ & ClusterCentroidBand' < 4)
%             cluster_probunits_window_thetaalpha(xi,evt_i) = sum(Sign_Taskid(evt_i,idxUnits == idx_ & ClusterCentroidBand' < 4),2)*100/SPCDensityUnits.nPairs(xi);
%             end
%             if sum(idxUnits == idx_ & ClusterCentroidBand'  == 4)
%             cluster_probunits_window_beta(xi,evt_i) = sum(Sign_Taskid(evt_i,idxUnits == idx_ & ClusterCentroidBand'  == 4),2)*100/SPCDensityUnits.nPairs(xi);
%             end
%         end
%     end


% sign_rel = nan(1,height(SPCDensityUnits));
% r_rel = nan(1,height(SPCDensityUnits));
% pperm_rel = nan(1,height(SPCDensityUnits));

% change IFR vs t-SPC change (beta desynch)

[r_rel,~,pperm_rel,fh] = corr_permute( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),BetaDesynch(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),BetaDesynch(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')

ylabel('beta desynch t-SPC density change (rest - speech)')
xlabel('IFR change (z-score w.r.t. rest) [spks/s]')
ylim([-11 53])
axfh.Parent=tcl;
axfh.Layout.Tile=2;

saveFigures(fh,fullfile(PATHS.saveFigures,'BetaDesynch_vs_IFR'));


% compare IFR change in units with change of SPC or not
IFR_noSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & BetaDesynch == 0);
IFR_yesSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & BetaDesynch ~= 0);

p = permutation_diffTest_indsamples(IFR_noSPC_change, IFR_yesSPC_change,500);

figure('position',[300 300 500 400])
tiledlayout(2,1)
nexttile
raincloud_plot(IFR_noSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_noSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_noSPC_change),iqr(IFR_noSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta desynch (no IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])
nexttile
raincloud_plot(IFR_yesSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_yesSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_yesSPC_change),iqr(IFR_yesSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta desynch (IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])


[r_rel,~,pperm_rel,fh] = corr_permute(SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),"Spearman",500,0);
%[r_rel,~,pperm_rel,fh] = corr_permute(SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0),BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0 & SPCDensityUnits.beta_speech == 0),"Spearman",500,0);

% change IFR vs t-SPC change (theta/alpha synch)

[r_rel,~,pperm_rel,fh] = corr_permute( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),BetaDesynch(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;
fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
%scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')
ylim([-3 64])
ylabel('thetaalpha synch t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

saveFigures(fh,fullfile(PATHS.saveFigures,'ThetaAlphaSynch_vs_IFR'));



% compare IFR change in units with change of SPC or not
IFR_noSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & ThetaAlphaSynch == 0);
IFR_yesSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & ThetaAlphaSynch ~= 0);
p = permutation_diffTest_indsamples(IFR_noSPC_change, IFR_yesSPC_change,500);

figure('position',[300 300 500 400])
tiledlayout(2,1)
nexttile
raincloud_plot(IFR_noSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_noSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_noSPC_change),iqr(IFR_noSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC theta/alpha desynch (no IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])
nexttile
raincloud_plot(IFR_yesSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_yesSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_yesSPC_change),iqr(IFR_yesSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC theta/alpha desynch (IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])

[r_rel,~,pperm_rel,fh] = corr_permute(SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),"Spearman",500,0);

% change IFR vs t-SPC change (beta rebound)

[r_rel,~,pperm_rel,fh] = corr_permute( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),BetaRebound(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean)),BetaRebound(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
scatter( SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaRebound(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')
ylim([-11 53])

ylabel('beta rebound t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

saveFigures(fh,fullfile(PATHS.saveFigures,'BetaRebound_vs_IFR'));


% compare IFR change in units with change of SPC or not
IFR_noSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & BetaRebound == 0);
IFR_yesSPC_change = SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & BetaRebound ~= 0);
p = permutation_diffTest_indsamples(IFR_noSPC_change, IFR_yesSPC_change,500);

figure('position',[300 300 500 400])
tiledlayout(2,1)
nexttile
raincloud_plot(IFR_noSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_noSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_noSPC_change),iqr(IFR_noSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta rebound (no IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])
nexttile
raincloud_plot(IFR_yesSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_yesSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_yesSPC_change),iqr(IFR_yesSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta rebound (IFR change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])


[r_rel,~,pperm_rel,fh] = corr_permute(SPCDensityUnits.FRzscore(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaRebound(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),"Spearman",500,0);


% change t-sPC beta desynch vs t-SPC beta rebound

[r_rel,~,pperm_rel,fh] = corr_permute( BetaDesynch(~isnan(SPCDensityUnits.FRmean)),BetaRebound(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( BetaDesynch(~isnan(SPCDensityUnits.FRmean)),BetaRebound(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
scatter( BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaRebound(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')
ylim([-11 53])
xlim([-11 53])

ylabel('beta rebound t-SPC density change (rest - speech)')
xlabel('beta desynch t-SPC density change (rest - speech)')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

saveFigures(fh,fullfile(PATHS.saveFigures,'BetaDesynch_vs_BetaRebound'));

[r_rel,~,pperm_rel,fh] = corr_permute(BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),BetaRebound(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),"Spearman",500,0);



% change t-sPC thetaalpha synch vs t-SPC beta desynch

[r_rel,~,pperm_rel,fh] = corr_permute( BetaDesynch(~isnan(SPCDensityUnits.FRmean)),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( BetaDesynch(~isnan(SPCDensityUnits.FRmean)),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
scatter( BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')
ylim([-3 64])
xlim([-11 53])

ylabel('thetaalpha t-SPC density change (rest - speech)')
xlabel('beta desynch t-SPC density change (rest - speech)')
axfh.Parent=tcl;
axfh.Layout.Tile=2;


% compare IFR change in units with change of SPC or not
IFR_noSPC_change = BetaRebound(~isnan(SPCDensityUnits.FRmean) & ThetaAlphaSynch == 0);
IFR_yesSPC_change = BetaRebound(~isnan(SPCDensityUnits.FRmean) & ThetaAlphaSynch ~= 0);

p = permutation_diffTest_indsamples(IFR_noSPC_change, IFR_yesSPC_change,500);

figure('position',[300 300 500 400])
tiledlayout(2,1)
nexttile
raincloud_plot(IFR_noSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_noSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_noSPC_change),iqr(IFR_noSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta desynch (no theta change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])
nexttile
raincloud_plot(IFR_yesSPC_change,"box_on",1,"color",[.8 .8 .8],"alpha",.3)
hold on
xline(median(IFR_yesSPC_change),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(IFR_yesSPC_change),iqr(IFR_yesSPC_change)),'LabelOrientation','horizontal','linewidth',1.5)
box off
grid off
xlabel(" SPC beta desynch (theta change) [s]" )
ylabel(" PDF [a.u.] ")
yticks([0 0.15])
xlim([-10 30])



saveFigures(fh,fullfile(PATHS.saveFigures,'ThetaAlphaSynch_vs_BetaDesynch'));

[r_rel,~,pperm_rel,fh] = corr_permute(BetaDesynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),"Spearman",500,0);

% change t-sPC thetaalpha synch vs t-SPC beta rebound

[r_rel,~,pperm_rel,fh] = corr_permute( BetaRebound(~isnan(SPCDensityUnits.FRmean)),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

fh = figure('Position',[200 200 800 400]);
tcl = tiledlayout(1,2);
nexttile
scatter( BetaRebound(~isnan(SPCDensityUnits.FRmean)),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean)),'k','filled')

hold on
scatter( BetaRebound(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),ThetaAlphaSynch(~isnan(SPCDensityUnits.FRmean) & SPCDensityUnits.beta_ITI > 0),'r','filled')
ylim([-3 64])
xlim([-11 53])

ylabel('thetaalpha t-SPC density change (rest - speech)')
xlabel('beta rebound t-SPC density change (rest - speech)')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

saveFigures(fh,fullfile(PATHS.saveFigures,'ThetaAlphaSynch_vs_BetaRebound'));


%%
% base IFR vs t-SPC base
[r_rel,~,pperm_rel] = corr_permute( cluster_probunits_window(~isnan(SPCDensityUnits.FRbase),1),SPCDensityUnits.FRbase(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);

% change IFR vs t-SPC change
[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window(~isnan(SPCDensityUnits.FRbase),4),SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),cluster_probunits_window(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window(~isnan(SPCDensityUnits.FRbase),4),'k','filled')
ylabel('t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

% change IFR vs t-SPC change (beta desynch)

[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),'k','filled')
ylabel('beta desynch t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

% change IFR vs t-SPC change (beta reb)
[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),5) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1),SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),5) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1),'k','filled')
ylabel('beta rebound t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

% change IFR vs t-SPC change (alpha)
[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),4) - cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),1),SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(SPCDensityUnits.FRchange(~isnan(SPCDensityUnits.FRmean)),cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),4) - cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),1),'k','filled')
ylabel('that/alpha sync t-SPC density change (rest - speech)')
xlabel('IFR change (% change w.r.t. rest) [spks/s]')
axfh.Parent=tcl;
axfh.Layout.Tile=2;


% change alpha t-sCP change vs beta t-SPC change 
[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),4) - cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),1),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),4) - cluster_probunits_window_thetaalpha(~isnan(SPCDensityUnits.FRbase),1),cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),'k','filled')
ylabel('beta desynch t-SPC density change (rest - speech)')
xlabel('that/alpha sync t-SPC density change (rest - speech)')
axfh.Parent=tcl;
axfh.Layout.Tile=2;


% change beta desynch t-sCP change vs beta t-SPC rebound 
[r_rel,~,pperm_rel,fh] = corr_permute( cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),5) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1),"Spearman",500,0);
axfh = fh.CurrentAxes;

figure('Position',[200 200 800 400])
tcl = tiledlayout(1,2);
nexttile
scatter(cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),4),cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),5) - cluster_probunits_window_beta(~isnan(SPCDensityUnits.FRbase),1),'k','filled')
ylabel('beta desynch t-SPC density change (rest - speech)')
xlabel('beta rebound t-SPC density change (rest - speech)')
axfh.Parent=tcl;
axfh.Layout.Tile=2;

%% check SPCdensity across unit types

for fi = 2 : 6
    tmp = SPCDensityUnits{:,fi+4};
    for ui = 1 : numel(UnitOfInt)
        fprintf('%s type: %1.2f (%1.2f) %s-SPC density \n', ...
            UnitLabel(ui),mean(tmp(SPCDensityUnits.FRmod == UnitOfInt(ui)),'omitnan'), ...
            std(tmp(SPCDensityUnits.FRmod == UnitOfInt(ui)),[],'omitnan'), bands.name{fi});
    end

end


for fi = 2:6
end

