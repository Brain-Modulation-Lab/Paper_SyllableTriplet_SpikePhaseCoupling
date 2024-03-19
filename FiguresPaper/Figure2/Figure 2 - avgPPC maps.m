clc
clear all
close all

%% Step 0: settings of the script
Figure_2_settings;

%% Step 1. gather all SPC (pairs & clusters) information
DB = get_SPCdata(SUBJECTS,'all', PATHS);

%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB.Pairs);
%% Step 3. Extract SPC maps 

% 3.1 (PPC z-score) (all pairs)

% 3.1a Show PPCz on Speech
T = [-2.5 2];
nPerms = 500;
cfg_fh.Baseline = [-Inf EvtTimes.Speech(2)];

SPCmap = get_SPCmap(DB.Pairs,'PPCz', T, ...
    Tres,cfg,'Speech','all','default',nPerms);

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-.02 .18];
cfg_fh.CMap = linspecer;

fh = plotter_SPCmap(SPCmap, cfg_fh)

if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-all_evt-speech');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-speech');
    saveFigures(fh, figname)
end

% 3.1b Show PPCz on Cue
T = [-0.5 4];
nPerms = 500;
cfg.Baseline = [-Inf EvtTimes.Cue(2)];
SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Cue','all','default',nPerms);

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Cue;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-.02 .18];
cfg_fh.CMap = linspecer;
fh = plotter_SPCmap(SPCmap, cfg_fh)
if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-cue');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-all_evt-cue');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-all_evt-cue');
    saveFigures(fh, figname)
end

%% 
% 3.2 (PPC z-score) (only sign pairs)
% 3.2a Show PPCz on Speech
T = [-2.5 2];
nPerms = 500;
cfg.Baseline = [-Inf EvtTimes.Speech(2)];

SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','OnlySign','default',nPerms);
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CMap = linspecer;
% cfg_fh.CLim = [-.02 .18];
cfg_fh.NClustplot = 4;

fh = plotter_SPCmap(SPCmap, cfg_fh)

% figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-speech');
% saveFigures(fh, figname)
%         
if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-speech');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-onlysign_evt-speech');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-speech');
    saveFigures(fh, figname)
end



% 3.2b Show PPCz on Cue
T = [-0.5 4];
cfg.Baseline = [-Inf EvtTimes.Cue(2)];
nPerms = 500;
SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Cue','OnlySign','default',nPerms);
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Cue;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CMap = linspecer;
fh = plotter_SPCmap(SPCmap, cfg_fh)

%figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-cue');
if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-cue');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-onlysign_evt-cue');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-onlysign_evt-cue');
    saveFigures(fh, figname)
end

%% top connectivity pairs per unit 
% 3.3a (PPC z-score) (Speech)
T = [-2.5 2];
cfg.Baseline = [-Inf EvtTimes.Speech(2)];
nPerms = 500;

SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','Top','default',nPerms);
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CMap = linspecer;
cfg_fh.NClustplot = 15;
%cfg_fh.CLim = [-.3 .6];

fh = plotter_SPCmap(SPCmap, cfg_fh)
% figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-speech');
% saveFigures(fh, figname)
if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-speech');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-top_evt-speech');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-speech');
    saveFigures(fh, figname)
end


% 3.3b (PPC z-score) (Cue)
T = [-0.5 4];
cfg.Baseline = [-Inf EvtTimes.Cue(2)];
nPerms = 500;

SPCmap = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Cue','Top','default',nPerms);
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Cue;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CMap = linspecer;
cfg_fh.CLim = [-.3 .6];

fh = plotter_SPCmap(SPCmap, cfg_fh)
% figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-cue');
% saveFigures(fh, figname)

if nPerms > 0 & numel(fh) == 2
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-cue');
    saveFigures(fh{1}, figname)
    figname = fullfile(PATHS.saveFigures,'PPCtstat_pairs-top_evt-cue');
    saveFigures(fh{2}, figname)
else
    figname = fullfile(PATHS.saveFigures,'PPCzmap_pairs-top_evt-cue');
    saveFigures(fh, figname)
end

