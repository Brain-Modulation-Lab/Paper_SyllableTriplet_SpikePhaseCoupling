%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
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

TYPE_DB = 'default';
% PATH_FIGURES = ;



fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

PATHS = struct();
PATHS.DataDir      = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results_new', TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures_new',TYPE_DB);


%% plot Unit6 - STG, DBS3008, theta coupling example
SubjectOfInt = 'DBS3008';
EROIOfInt = 'G_temp_sup-Lateral L';
UnitOfInt = 'Unit6';
freqsOfInt = [8 12];


subj_i = find(contains(SUBJECTS,SubjectOfInt));

SUBJECT = SUBJECTS{subj_i};
fprintf('Starting computing PLV results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
PATH_E = fullfile(PATHS.DataDir ,SUBJECT,'E.mat');
PATH_S = fullfile(PATHS.DataDir ,SUBJECT,'S.mat');

PATH_MNI = fullfile(PATHS.DataDir,SUBJECT,'MNI.mat');
%PATH_STA = fullfile(PATH_RESULTS,SUBJECT,'STA.mat');

PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');


% load E
load(PATH_E)

% load S
load(PATH_S)

% load MNI
load(PATH_MNI)

% compute PLV;
%%


WIN_GET_SPIKENR = cfg.plv.WIN_GET_SPIKENR; % [sec], if you do not have lots of trials and thus not so many spikes you will want to make this bigger
tWidth_avgSpike = cfg.plv.tWidth_avgSpike;  % [sec] window size and..
tOffset_avgSpike = cfg.plv.tOffset_avgSpike;  % [sec] ...offset used to calculate the average spike number

NUM_PERMS = cfg.plv.NUM_PERMS;       % number of permutations, set at least to 200
PERM_WIN_LEN_SEC = cfg.plv.PERM_WIN_LEN_SEC; % this is used to permute the phase-signal, the exact number does not matter,
% it only needs to be longer than the largest bin that you'll use to compute the PLV with


LENGTH_WINDOW = cfg.plv.LENGTH_WINDOW;  % this is to set NBINS_WINDOW of bins in a window LENGTH_WINDOW
NBINS_WINDOW = cfg.plv.NBINS_WINDOW;
THR_NAN = cfg.plv.THR_NAN;
MIN_TRIALS = cfg.plv.MIN_TRIALS;
MAX_NAN = cfg.plv.MAX_NAN;

freqsOfInt = cfg.plv.freqsOfInt;   % Define your frequencies of interest

% find channel labels
E_channel = E.label(contains(MNI.E.atlas_label_Destrieux,EROIOfInt));
S_channel = UnitOfInt;

%
% select neuron
idx_neuron = find(strcmpi(S.label,S_channel));
spikeTimes = S.timestamp{idx_neuron};
% need coloumn vector for spike times
if size(spikeTimes,1) == 1
    spikeTimes = spikeTimes';
end

% spike times in seconds
try
    sess_neuron = unique(S.epochs{idx_neuron}.session_id_coding);
catch
    sess_neuron = unique(S.epochs{idx_neuron}.session_id);
end
% just a check
assert(numel(sess_neuron) == 1, "Single-Unit is present during more ECoG sessions - need to handtune the code!")


% select E channel  during proper session
ch_i = 2;
cfg_ft = [];
cfg_ft.channel = {E_channel{ch_i}};
E_data = ft_selectdata(cfg_ft,E);
E_data = E_data.trial{sess_neuron};

nan_perc = 1E2*sum(isnan(E_data))/numel(E_data);
% select sampling rate
SR = E.fsample;

% get neurons events (need to correct for non-nan ecog epochs, it could be brutal the correction losing so many trials...)
epochs = S.epochs{idx_neuron};
% find epochs with nan in E_data
time_nan = E.time{sess_neuron}(isnan(E_data));
try
    flag_nan = sum(time_nan >= epochs.starts_coding & time_nan <= epochs.ends_coding,2) > THR_NAN;
catch
    flag_nan = sum(time_nan >= epochs.starts & time_nan <= epochs.ends,2) > THR_NAN;
end
epochs_valid = epochs(~flag_nan,:);

% if spikes lands in
ntrials = height(epochs_valid);



% adjust spike-times and events by first E timing
offTime = E.time{sess_neuron}(1);
spikeTimes = spikeTimes - offTime;

Evt = table();
Evt.cue_onset = epochs_valid.stim1_onset - offTime;
Evt.cue_offset = epochs_valid.stim3_offset - offTime;
Evt.speech_onset = epochs_valid.syl1_onset - offTime;
Evt.speech_offset = epochs_valid.syl3_offset - offTime;




%     currFile = fileNames{f};
%     load([Paths.DataDir, currFile,'.mat'])
%     Evt = data.Evt;
%     LFP = data.LFP;
%     SR  = data.LFP_SR;
%     spikeTimes = data.spikeTimes;  % spike times in seconds

%% Load spikes
spike_Vect = zeros(1, length(E_data));
idcs_tmp = round(spikeTimes*SR);
idcs_tmp(idcs_tmp == 0) = []; % remove the index of 0 if there was one
spike_Vect(idcs_tmp) = 1;

%% Get epochs around your main event of interest to calculate the number of spikes you want to include per window
% based on the average within a window of a duration defined in WIN_GET_SPIKENR
[spike_atSpeech_trials]  =  getTrialsMatrix(Evt.cue_onset, spike_Vect, 1/SR, tWidth_avgSpike, tOffset_avgSpike);  % get epochs, aligned to your main event of interest
spikeTrialSum = sum(spike_atSpeech_trials,2);
nSpikes       = round(sum(spikeTrialSum) / tWidth_avgSpike * WIN_GET_SPIKENR);


%% Compute the time points of interest for each trial
allEvts = [];
trls_speech = ~isnan(Evt.cue_onset + Evt.cue_offset + Evt.speech_onset + Evt.speech_offset);
%trls_mvmt(end) = false; % remove the last event because the data was too short in some cases, @Matteo, you can delete this line

allEvts(1,:) = Evt.cue_onset(trls_speech) - 0.75;
allEvts(2,:) = Evt.cue_onset(trls_speech);
allEvts(3,:) = Evt.cue_offset(trls_speech);
allEvts(4,:) = Evt.speech_onset(trls_speech);
allEvts(5,:) = Evt.speech_offset(trls_speech);
allEvts(6,:) = Evt.speech_offset(trls_speech) + 0.75;

num_bins = round(diff(mean(allEvts,2)) / LENGTH_WINDOW * NBINS_WINDOW); % decide how many bins you want to have, this here is set to create 10 bins per 0.5 seconds

winCenters = [];
for e = 1:(size(allEvts,1)-1)
    winCenters_tmp = [];
    for trl = 1:size(allEvts,2)
        % Get the appropriate centers for each trial by
        % making equidistant points between neighbouring
        % events
        winCenters_tmp(trl,:) = linspace(allEvts(e,trl), allEvts(e+1,trl), num_bins(e));
    end
    if e == 1
        winCenters = [winCenters, winCenters_tmp]; % simply append the first time
    else % and for any subsequent events, remove the first winCenter point, as this is
        % identical to the last winCenter point of the previous event
        winCenters = [winCenters, winCenters_tmp(:, 2:end)];
    end
end


if cfg.flip.flip_check
    flip_flag = flip_polarity(E_data, SR, cfg.flip);

end

if cfg.plv.gamma_env
    [bFOI, aFOI] = butter(cfg.flip.ordFfoi, cfg.plv.gamma_env_BPF/(SR/2));
    tmp= E_data;
    tmp(isnan(E_data)) = 0;
    GammaLFP_filt = filtfilt(bFOI, aFOI,tmp) ;
    GammaLFP_filt = abs(hilbert(GammaLFP_filt));
end

freqOfInt = freqsOfInt;

%% Extract the LFP phase for the frequency of interest
%         LFP_filt1 = ft_preproc_highpassfilter(LFP, SR, freqOfInt(1), 4, 'but','twopass'); % 3rd order butterworth high-pass filter
%         LFP_filt1 = ft_preproc_lowpassfilter(LFP_filt1, SR, freqOfInt(2), 4, 'but','twopass'); % 3rd order butterworth low-pass filter
%         d = designfilt('bandpassiir', 'FilterOrder',4, ...
%             'HalfPowerFrequency1',freqOfInt(1),'HalfPowerFrequency2',freqOfInt(2), ...
%             'SampleRate',SR);
%         LFP_filt = filter(d, LFP);  %% !!! change this into a two-pass filter (passed forwards and backwards, so that you'll get zero phase distortion)
%         LFP_filt = filtfilt(d, LFP);  % watch out, this gives weird
%         results if you have the fieldtrip toolbox or EEGLAB toolbox added
%         to your path, currently not working

if ~cfg.plv.gamma_env
    [bFOI, aFOI] = butter(cfg.flip.ordFfoi, freqOfInt/(SR/2));
    tmp= E_data;
    tmp(isnan(E_data)) = 0;
    LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
    phase = angle(hilbert(LFP_filt));
else
    [bFOI, aFOI] = butter(cfg.flip.ordFfoi, freqOfInt/(SR/2));
    tmp = GammaLFP_filt;
    tmp(isnan(GammaLFP_filt)) = 0;
    LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
    phase = angle(hilbert(LFP_filt));
end

if flip_flag
    phase = calc_shiftPhase(phase,pi);
end
%phase(isnan(E_data)) = nan;
% check phase correct


%% Calculate PLV in the windows and estimate the binLength (or window width) for each window
idcs = round(winCenters*SR);
w = 59;
% Important function that finds the appropriate
% window-width to match the number of spikes across
% windows that are used to calculate the PLV
[curr_spks, store] = calc_binLen_selectSpikes(idcs, w, spike_Vect, nSpikes, store, nan, nan);
spks_inphase = find(curr_spks == 1);
stopspks = spks_inphase(581) : spks_inphase(621);
currPhases = phase(stopspks);
currLFP = LFP_filt(stopspks);
spks_inphase = spks_inphase(581:621);
spon = epochs_valid.syl1_onset;
spoff = epochs_valid.syl3_offset;
%%
for ii = 1:60%52 %14
    figure('renderer','painters','Position',[300 100 1200 400])
    idx = E.time{1} <=  spoff(ii)+1  & E.time{1} >=  spon(ii)-1;
    spks_toplot = S.timestamp{idx_neuron}(S.timestamp{idx_neuron} <= spoff(ii)+1  & S.timestamp{idx_neuron} >=  spon(ii)-1);
    idxspike = dsearchn(E.time{1}(idx)',spks_toplot');
    nexttile
    hold on
    plot([spon(ii) spon(ii)], [-6 6], 'r--')
    plot([spoff(ii) spoff(ii)], [-6 6], 'r--')
    plot(E.time{1}(idx),zscore(LFP_filt(idx)),'k')
    plot(repmat(spks_toplot,2,1),[(3.2)*ones(1,numel(spks_toplot)) ; (3.6)*ones(1,numel(spks_toplot))] ,'k')
    ylim([-6 6])
    yyaxis right
    ax = gca;
    ax.YColor = 'b';
    plot(E.time{1}(idx),phase(idx),'b')
    hold on
    tmp = phase(idx);
    scatter(spks_toplot,tmp(idxspike),20,'b','filled'   )
    ylim([-4 20])
    box off
    axis off
    pause()
end



% #########################
%     %% Save the file
%     if ~exist([Paths.saveMatFiles, '/PLV'])
%         mkdir([Paths.saveMatFiles, '/PLV'])
%     end
%     save([Paths.saveMatFiles, '\PLV\PLV_', currFile, '_', LFP_channel, '_', spike_channel, 'NUM_PERM=', num2str(NUM_PERMS), '_WIN=', sprintf('%.2f', WIN_GET_SPIKENR), '.mat'],  'store', '-v7.3')


figname =fullfile(PATH_PROJECT,'Supplementary_Analysis','Figures','SPC_alpha_snippet');
saveas(gcf,figname,'fig')
saveas(gcf,figname,'png')
saveas(gcf,figname,'eps')

