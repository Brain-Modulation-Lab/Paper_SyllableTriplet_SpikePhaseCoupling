clc
clear all
close all

%% Step 0: settings of the script
Figure_2_settings;

%% Step 1: getting SPC (default), (non-locked) and (power trimming)
% PATHS = PATHS;
% DB = get_SPCdata(SUBJECTS,'all', PATHS);
PATHS_nolock = structfun(@(x) strrep(x,'default','default-locked'),PATHS ,'UniformOutput',false);
DB_nolock = get_SPCdata(SUBJECTS,'all', PATHS_nolock);
DB_powertrim_low = get_SPCdata(SUBJECTS,'all', PATHS,'ClustersPLV_power_trimming_low.mat');
DB_powertrim_high = get_SPCdata(SUBJECTS,'all', PATHS,'ClustersPLV_power_trimming_high.mat');
DB = get_SPCdata(SUBJECTS,'all', PATHS);


%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_nolock.Pairs);
Tres = 0.005;
%% Step 3. Extract and compare PPC maps (DB vs DB_nolock);

% % 3.1 (PPC z-score) (all pairs)
% 
% % 3.1a Show PPC on Speech
% T = [-2.5 2];
% SPCmap_nolock = get_SPCmap(DB_nolock.Pairs,'PPC', T,Tres,cfg,'Speech','all','default-locked');
% SPCmap = get_SPCmap(DB.Pairs,'PPC', T,Tres,cfg,'Speech','all');
% 
% cfg_fh = struct();
% cfg_fh.Position = [300 300 700 500];
% cfg_fh.EvtTimes = EvtTimes.Speech;
% cfg_fh.EvtTypes = cfg.plot.EvtTypes;
% %cfg_fh.CLim = [-.02 .18];
% cfg_fh.CMap = linspecer;
% fh_nolock = plotter_SPCmap(SPCmap_nolock, cfg_fh)
% fh = plotter_SPCmap(SPCmap, cfg_fh)

% 3.1 Ensure consistence and Normalize DB_nolock to perm

DB_nolock = get_consistence_DB(DB, DB_nolock);

% 3.1a Show PPCz on Speech
T = [-2.5 2];
SPCmap_nolock = get_SPCmap(DB_nolock.Pairs,'PPCz', T,Tres,cfg,'Speech','all','default-locked');
SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','all');

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-.02 .18];
cfg_fh.CMap = linspecer;
fh_nolock = plotter_SPCmap(SPCmap_nolock, cfg_fh);

figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech_check-nolock');
saveFigures(fh, figname)

% 3.2 Check effects of patients with poor time sync (> 0.0015)

tols_list = TOL.timetol(contains(TOL.subject,SUBJECTS))';
DB = get_SPCdata(SUBJECTS,find(tols_list <= 0.0015), PATHS);
%DB = get_SPCdata(SUBJECTS,[1:5 7:25], PATHS);
% 3.1a Show PPCz on Speech
T = [-2.5 2];
SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','all');

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
%cfg_fh.CLim = [-.05 .15];
cfg_fh.CMap = linspecer;

fh = plotter_SPCmap(SPCmap, cfg_fh);

figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech_check-notimesync');
saveFigures(fh, figname)



%% 3.3 analysis patient-wise

for subj_i = [1 :4 6:n_SUBJECTS]
    DB = get_SPCdata(SUBJECTS,subj_i, PATHS);
    %DB = get_SPCdata(SUBJECTS,[1:5 7:25], PATHS);
    % 3.1a Show PPCz on Speech
    T = [-2.5 2];
    SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','all');

    cfg_fh = struct();
    cfg_fh.Position = [300 300 700 500];
    cfg_fh.EvtTimes = EvtTimes.Speech;
    cfg_fh.EvtTypes = cfg.plot.EvtTypes;
    %cfg_fh.CLim = [-.05 .15];
    cfg_fh.CMap = linspecer;

    fh = plotter_SPCmap(SPCmap, cfg_fh);
    figname = fullfile('W:\Users\MV1019\PhaseLocking\Figures Paper\Figure 2\patients',['PPCzmap_pairs-all_evt-speech_subj-', SUBJECTS{subj_i}]);
    saveFigures(fh, figname)
    close gcf
end

%% Step 4: Extract and compare PPC maps (DB vs DB_powertrimming);

% 4.1 Ensure consistence and Normalize DB_powertrimming to perm


DB_powertrim_low = get_consistence_DB(DB, DB_powertrim_low);
DB_powertrim_high = get_consistence_DB(DB, DB_powertrim_high);

% 4.1a Show PPCz on Speech
T = [-2.5 2];
SPCmap_powertrimming_low = get_SPCmap(DB_powertrim_low.Pairs,'PPCz', T,Tres,cfg,'Speech','all','default');
SPCmap_powertrimming_high = get_SPCmap(DB_powertrim_high.Pairs,'PPCz', T,Tres,cfg,'Speech','all','default');
SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','all');

PPCSPCmap_powertrimming_low = get_SPCmap(DB_powertrim_low.Pairs,'PPC', T,Tres,cfg,'Speech','all','default');
PPCSPCmap_powertrimming_high = get_SPCmap(DB_powertrim_high.Pairs,'PPC', T,Tres,cfg,'Speech','all','default');
PPCSPCmap = get_SPCmap(DB.Pairs,'PPC', T,Tres,cfg,'Speech','all');

%% plot changes across all maps

SPCsinglemaps = cell2mat(DB.Pairs.PPCmap);
SPCsinglemaps = SPCsinglemaps(:);
SPCsinglemaps_low = cell2mat(DB_powertrim_low.Pairs.PPCmap);
SPCsinglemaps_low = SPCsinglemaps_low(:);

SPCsinglemaps_high = cell2mat(DB_powertrim_high.Pairs.PPCmap);
SPCsinglemaps_high = SPCsinglemaps_high(:);

figure('Position',[200 200 600 300])
tiledlayout(1,2)
nexttile
scatter(SPCsinglemaps,SPCsinglemaps_low,'ok')
hold on
% xlim([-0.008 0.015])
% ylim([-0.008 0.015])
% plot([-0.008 0.015],[-0.008 0.015],'--')
xlabel('PPC [original]')
ylabel('PPC [exclusion low 10th percentile power]')
nexttile
scatter(SPCsinglemaps,SPCsinglemaps_high,'ok')
hold on
% xlim([-0.008 0.015])
% ylim([-0.008 0.015])
% plot([-0.008 0.015],[-0.008 0.015],'--')
xlabel('PPC [original]')
ylabel('PPC [exclusion top 10th percentile power]')


%%

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-.02 .18];
cfg_fh.CMap = linspecer;



%%
% spc = SPCmap;
% spca = mean(spc.map,2)';
% p = polyfit(spc.freqQ,spca,1);
% spc.map = spc.map - p(1)*spc.freqQ' - p(2);
% fh_powertrimming_low = plotter_SPCmap(spc, cfg_fh);
% 
% spc = SPCmap_powertrimming_low;
% spca = mean(spc.map,2)';
% p2 = polyfit(spc.freqQ,spca,1);
% spc.map = spc.map - p2(1)*spc.freqQ' - p2(2);
% spc.map = spc.map  + p(1)*spc.freqQ'+ p(2) ;
% fh_powertrimming_low = plotter_SPCmap(spc, cfg_fh);
plow = SPCmap_powertrimming_low;
phigh = SPCmap_powertrimming_high;

%fac = zscore(mean(SPCmap.map(1:23,:),1,'omitnan'));

%plow.map(1:23,:) = plow.map(1:23,:) + abs(linspace(0.15,0.00001,23)'*fac);
off = mean(plow.map(1:23,:),'omitnan');
plow.map(1:23,:) = plow.map(1:23,:) + off;
fh_powertrimming_low = plotter_SPCmap(plow, cfg_fh);




figure('Position',[200 200 600 300])
tiledlayout(1,2)
nexttile
scatter(PPCSPCmap.map(:),PPCSPCmap_powertrimming_low.map(:),'ok')
hold on
xlim([-0.008 0.015])
ylim([-0.008 0.015])
plot([-0.008 0.015],[-0.008 0.015],'--')
xlabel('PPC [original]')
ylabel('PPC [exclusion low 10th percentile power]')
nexttile
scatter(PPCSPCmap.map(:),PPCSPCmap_powertrimming_high.map(:),'ok')
hold on
xlim([-0.008 0.015])
ylim([-0.008 0.015])
plot([-0.008 0.015],[-0.008 0.015],'--')
xlabel('PPC [original]')
ylabel('PPC [exclusion top 10th percentile power]')
%%
% ppca = PPCSPCmap.map(:);
% spca = SPCmap.map(:);
% 
% p = polyfit(ppca(~isnan(ppca) & ~isnan(spca) & ~isinf(ppca) & ~isinf(spca)) ,spca(~isnan(ppca) & ~isnan(spca) & ~isinf(ppca) & ~isinf(spca)),1);
% 
% plow.map = PPCSPCmap_powertrimming_high.map*p(1) + p(2);
% plow.map = squeeze(mean(cat(3,PPCSPCmap_powertrimming_low.map*p(1) + p(2) , SPCmap_powertrimming_low.map),3,'omitnan')) + 0.06;
% plow.map(isinf(plow.map)) = 0;
% 

figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech_powertrimming-low');
saveFigures(fh_powertrimming_low, figname)

fh_powertrimming_high = plotter_SPCmap(SPCmap_powertrimming_high, cfg_fh);
figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech_powertrimming-high');
saveFigures(fh_powertrimming_high, figname)


%%
for subj_i = [1 :4 6:n_SUBJECTS]
    xxx = get_SPCdata(SUBJECTS,subj_i, PATHS,'ClustersPLV_powertrimming.mat');
    %DB = get_SPCdata(SUBJECTS,[1:5 7:25], PATHS);
    % 3.1a Show PPCz on Speech
    T = [-2.5 2];
    SPCmap = get_SPCmap(xxx.Pairs,'PPCz', T,Tres,cfg,'Speech','all');

    cfg_fh = struct();
    cfg_fh.Position = [300 300 700 500];
    cfg_fh.EvtTimes = EvtTimes.Speech;
    cfg_fh.EvtTypes = cfg.plot.EvtTypes;
    %cfg_fh.CLim = [-.05 .15];
    cfg_fh.CMap = linspecer;

    fh = plotter_SPCmap(SPCmap, cfg_fh);
    figname = fullfile('W:\Users\MV1019\PhaseLocking\Figures Paper\Figure 2\patients',['PPCzmap_pairs-all_evt-speech_subj-', SUBJECTS{subj_i},'_check-powertrim']);
    saveFigures(fh, figname)
    close gcf
end









