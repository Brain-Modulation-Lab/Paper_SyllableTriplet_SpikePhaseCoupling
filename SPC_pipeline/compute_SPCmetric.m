function SPC = compute_SPCmetric(E, S, cfg_analysis)

% INPUT:
%           E (LFP channel): struct with fields:
%                                       data : 1 x N (N = # samples)
%                                       time : 1 x N (N = # samples)
%                                       epoch: Ev x 3 (Ev = # events, onset, left offset, right offset)
%                                              final events: [epoch(:,1) - epoch(:,2) epoch(:,1) + epoch(:,3)]
%           S (spike timestamps in [s]): 1 X M (timestamps of spikes, M = #
%           spikes)
%           cfg_analysis: struct with fields (output of set_configs.m)
% OUTPUT:
%           SPC: struct that contaons the spike-phase coupling metrics
%

% timestamps must be corrected with first time of LFP;
WIN_GET_SPIKENR = cfg_analysis.plv.WIN_GET_SPIKENR; % [sec], if you do not have lots of trials and thus not so many spikes you will want to make this bigger
tWidth_avgSpike = cfg_analysis.plv.tWidth_avgSpike;  % [sec] window size and..
tOffset_avgSpike = cfg_analysis.plv.tOffset_avgSpike;  % [sec] ...offset used to calculate the average spike number

NUM_PERMS = cfg_analysis.plv.NUM_PERMS;       % number of permutations, set at least to 200
PERM_WIN_LEN_SEC = cfg_analysis.plv.PERM_WIN_LEN_SEC; % this is used to permute the phase-signal, the exact number does not matter,
% it only needs to be longer than the largest bin that you'll use to compute the PLV with


LENGTH_WINDOW = cfg_analysis.plv.LENGTH_WINDOW;  % this is to set NBINS_WINDOW of bins in a window LENGTH_WINDOW
NBINS_WINDOW = cfg_analysis.plv.NBINS_WINDOW;
THR_NAN = cfg_analysis.plv.THR_NAN;
MIN_TRIALS = cfg_analysis.plv.MIN_TRIALS;
MAX_NAN = cfg_analysis.plv.MAX_NAN;

freqsOfInt = cfg_analysis.plv.freqsOfInt;   % Define your frequencies of interest


store = struct();

% LFP signal
E_data = E.data;
E_time = E.time;
% Neuron spiketimes
spikeTimes = S;
% Events
Events = E.epoch(:,1) + [- E.epoch(:,2) E.epoch(:,3)];

% need coloumn vector for spike times
if size(spikeTimes,1) == 1
    spikeTimes = spikeTimes';
end

% correct if locked
if cfg_analysis.locked
    x_ = [];
    for ee = 1 : size(Events,1)
        start_ = Events(1, ee);
        stop_ = Events(2,ee);
        idx_ = E_time >= start_ & E_time <= stop_;
        x_ = [x_ ;E_data(1, idx_)];
    end
    evoked = mean(x_,'omitnan');
    for ee = 1 : size(Events,1)
        start_ = Events(1, ee);
        stop_ = Events(2,ee);
        idx_ = E_time >= start_ & E_time <= stop_;
        E_data(idx_) = E_data(idx_) - evoked;
    end
end


% adjust spike-times and events by first E timing
offTime = E_time(1);
spikeTimes = spikeTimes - offTime;


% modify from here for publishing the code to the repository
% ------------------------------------------
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

if cfg_analysis.flip.flip_check
    flip_flag = flip_polarity(E_data, SR, cfg_analysis.flip);

end

if cfg_analysis.plv.gamma_env
    [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, cfg_analysis.plv.gamma_env_BPF/(SR/2));
    tmp= E_data;
    tmp(isnan(E_data)) = 0;
    GammaLFP_filt = filtfilt(bFOI, aFOI,tmp) ;
    GammaLFP_filt = abs(hilbert(GammaLFP_filt));
end


for frq = 1:size(freqsOfInt,1)
    freqOfInt = freqsOfInt(frq,:);

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

    if ~cfg_analysis.plv.gamma_env
        [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, freqOfInt/(SR/2));
        tmp= E_data;
        tmp(isnan(E_data)) = 0;
        LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
        phase = angle(hilbert(LFP_filt));
    else
        [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, freqOfInt/(SR/2));
        tmp = GammaLFP_filt;
        tmp(isnan(GammaLFP_filt)) = 0;
        LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
        phase = angle(hilbert(LFP_filt));
    end

    % handle power trim
    if  cfg_analysis.powertrim
        envtmp = abs(hilbert(LFP_filt));
        switch  cfg_analysis.powertrim_direction
            case 'low'
                envtmp_low = envtmp <= prctile(envtmp(~isnan(E_data)),cfg_analysis.powertrim_prctile);
            case 'high'
                envtmp_low = envtmp >= prctile(envtmp(~isnan(E_data)),cfg_analysis.powertrim_prctile);
        end
        envtmp_low = setdiff(find(envtmp_low),find(isnan(E_data)));
        spike_Vect_foi = spike_Vect;
        spike_Vect_foi(envtmp_low) = 0;
    else
        spike_Vect_foi = spike_Vect;
    end

    if flip_flag
        phase = calc_shiftPhase(phase,pi);
    end
    %phase(isnan(E_data)) = nan;
    % check phase correct


    %% Calculate PLV in the windows and estimate the binLength (or window width) for each window
    idcs = round(winCenters*SR);
    fprintf('Pair %d: %s-%s @freq: %i\n',pair_i, E_channel, S_channel, freqOfInt(1)) %currFile, LFP_channel, spike_channel, freqOfInt(1));
    for w = 1:size(winCenters,2)
        % Important function that finds the appropriate
        % window-width to match the number of spikes across
        % windows that are used to calculate the PLV
        [curr_spks, store] = calc_binLen_selectSpikes(idcs, w, spike_Vect_foi, nSpikes, store, nan, nan);

        currPhases                = phase(curr_spks==1);
        [PLV_original, meanAngle] = comp_PLV(currPhases);
        store.PLV(frq, w)         = PLV_original;
        store.phase(frq, w)       = meanAngle;
        store.PPC(frq, w)          = correct_PLV(PLV_original, sum(~isnan(currPhases))); % correct plv
        store.ES(frq, w)           = calc_effectsizePPC(store.PPC(frq, w)); % correct effect size

        % Get the phase for the permutation procedure:
        %phase(isnan(phase)) = 0;
        [phases_atSpeech_trls, phase_atSpeech_vect]  = getTrialsMatrix(winCenters(:,w), phase, 1/SR,PERM_WIN_LEN_SEC, PERM_WIN_LEN_SEC/2);
        for perm = 1:NUM_PERMS
            rng(perm);  % use rng here, otherwise the perm distribution
            % would be re-shuffled for each window/frequency
            randIdcs = randperm(size(phases_atSpeech_trls,2));  % create random indices to shuffle the order of trials
            shuffledPhaseTrials = phases_atSpeech_trls(:,randIdcs);
            tmp              = phase_atSpeech_vect;       % this is the full time series, but only includes
            tmp(~isnan(tmp)) = shuffledPhaseTrials(:);  % refill the long vector with the shuffled trials
            shuffledPhase    = tmp;
            [PLV_perm, meanAngle_perm]     = comp_PLV(shuffledPhase(curr_spks==1));
            store.PLV_perm(frq, w, perm)   = PLV_perm;
            store.phase_perm(frq, w, perm) = meanAngle_perm;
            store.PPC_perm(frq, w, perm)   = correct_PLV(PLV_perm, sum(~isnan(shuffledPhase(curr_spks==1))));
        end
    end
end





