%% this script checks for beta and gamma modulation

clc
clear all
close all

PATH_SERVER = 'Z:\DBS';
PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
PATH_OUTPUT = fullfile(PATH_PROJECT,'Supplementary_Analysis');

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

fmin  = 5;
fmax = 180;
nfreq = 70;

fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

TYPE_DB = 'default';
PATHS = struct();
PATHS.DataDir      = fullfile('W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results');
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


%%

events = {'stim1_onset','syl1_onset'};
toDo = 6; %[1:4 6:n_SUBJECTS];
for subj_i = toDo

    % TF = [];
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Starting computing TF results for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    PATH_E = fullfile(PATHS.DataDir ,SUBJECT,'E.mat');

    %     PATH_MNI = fullfile(PATH_RESULTS,SUBJECT,'MNI.mat');
    %PATH_STA = fullfile(PATH_RESULTS,SUBJECT,'STA.mat');

    PATH_ANNOT = fullfile(PATH_SERVER,SUBJECT,'Preprocessed Data','Sync','annot');


    % load E
    load(PATH_E)
    E_channel = E.label;
    nEChannels = numel(E_channel);
    try
        sessions = unique(E.epochs.session_id_coding)';
    catch
        sessions = unique(E.epochs.session_id)';
    end

    [E_cuelocked, E_cuelocked_avg, E_cuelocked_detr, ...
        E_speechlocked, E_speechlocked_avg, E_speechlocked_detr] = deal(cell(1,max(sessions)));
    for sess_i = sessions
        try
            epochs = E.epochs(E.epochs.session_id_coding == sess_i,:);
        catch
            epochs = E.epochs(E.epochs.session_id == sess_i,:);
        end

        for event_i = 1 : 2
            switch event_i
                case 1% cue
                    cfg = [];
                    cfg.epoch = epochs;

                    %         cfg.timelock = evts_t0{1};
                    %         E_cuelocked{sess_i} = bml_redefinetrial(cfg,E);
                    cfg.timelock = events{event_i};
                    E_cuelocked{sess_i} = bml_redefinetrial(cfg,E);
                    % calculate trial average
                    %         E_cuelocked_avg{sess_i} = E_cuelocked{sess_i};
                    %         E_cuelocked_avg{sess_i}.time = E_cuelocked_avg{sess_i}.time{1};
                    %         E_cuelocked_avg{sess_i}.trial = {squeeze(mean(cat(3,E_cuelocked{sess_i}.trial{:}),3,'omitnan'))};
                    E_cuelocked_avg{sess_i} = E_cuelocked{sess_i};
                    E_cuelocked_avg{sess_i}.time = E_cuelocked_avg{sess_i}.time{1};
                    E_cuelocked_avg{sess_i}.trial = {squeeze(mean(cat(3,E_cuelocked{sess_i}.trial{:}),3,'omitnan'))};
                    % detrend trial average
                    %         E_cuelocked_detr{sess_i} =   E_cuelocked{sess_i};
                    %         E_cuelocked_detr{sess_i}.trial = cellfun(@(x) x - E_cuelocked_avg{sess_i}.trial{1,1}, E_cuelocked{sess_i}.trial,'uni',0);
                    E_cuelocked_detr{sess_i} =   E_cuelocked{sess_i};
                    E_cuelocked_detr{sess_i}.time = E_cuelocked_detr{sess_i}.time{1};
                    E_cuelocked_detr{sess_i}.trial = cellfun(@(x) x - E_cuelocked_avg{sess_i}.trial{1,1}, E_cuelocked{sess_i}.trial,'uni',0);

                    % calculate TF
                    disp(' TF decomposition in locked activity' )
                    cfg = [];
                    cfg.times = E_cuelocked_avg{sess_i}.time;
                    cfg.trials = 1;
                    cfg.pnts = numel(E_cuelocked_avg{sess_i}.time);
                    cfg.srate = E_cuelocked_avg{sess_i}.fsample;
                    cfg.trials = numel(E_cuelocked{sess_i}.trial);
                    E_cuelocked_avg{sess_i}.TF = wav_decomp(E_cuelocked_avg{sess_i}.trial{1,1},cfg,fmin,fmax,nfreq,'means',E_cuelocked_avg{sess_i}.time);


                    disp(' TF decomposition in (locked + non-locked) activity' )
                    cfg = [];
                    cfg.times = E_cuelocked{sess_i}.time{1,1};
                    cfg.trials = 1;
                    cfg.pnts = numel(E_cuelocked{sess_i}.time{1,1});
                    cfg.srate = E_cuelocked{sess_i}.fsample;
                    cfg.trials = numel(E_cuelocked{sess_i}.trial);
                    E_cuelocked{sess_i}.TF = wav_decomp(cat(3,E_cuelocked{sess_i}.trial{:}),cfg,fmin,fmax,nfreq,'means',E_cuelocked{sess_i}.time{1,1});

                   % calculate TF baseline norm
                    baselinetime = [(epochs.stim1_onset - .6) - epochs.stim1_onset  (epochs.stim1_onset - .1) - epochs.stim1_onset];
                    E_cuelocked{sess_i}.TFnorm = wav_decomp(E_cuelocked{sess_i}.trial{1,1},cfg,fmin,fmax,nfreq,'means',E_cuelocked{sess_i}.time{1,1},baselinetime);



                    
                    disp(' TF decomposition in non-locked activity' )
                    cfg = [];
                    cfg.times = E_cuelocked_detr{sess_i}.time;
                    cfg.trials = 1;
                    cfg.pnts = numel(E_cuelocked_detr{sess_i}.time);
                    cfg.srate = E_cuelocked_detr{sess_i}.fsample;
                    cfg.trials = numel(E_cuelocked_detr{sess_i}.trial);
                    E_cuelocked_detr{sess_i}.TF = wav_decomp(cat(3,E_cuelocked_detr{sess_i}.trial{:}),cfg,fmin,fmax,nfreq,'means',E_cuelocked_detr{sess_i}.time);


                case 2 % speech
                    cfg = [];
                    cfg.epoch = epochs;

                    %         cfg.timelock = evts_t0{1};
                    %         E_cuelocked{sess_i} = bml_redefinetrial(cfg,E);
                    cfg.timelock = events{event_i};
                    E_speechlocked{sess_i} = bml_redefinetrial(cfg,E);
                    % calculate trial average
                    %         E_cuelocked_avg{sess_i} = E_cuelocked{sess_i};
                    %         E_cuelocked_avg{sess_i}.time = E_cuelocked_avg{sess_i}.time{1};
                    %         E_cuelocked_avg{sess_i}.trial = {squeeze(mean(cat(3,E_cuelocked{sess_i}.trial{:}),3,'omitnan'))};
                    E_speechlocked_avg{sess_i} = E_speechlocked{sess_i};
                    E_speechlocked_avg{sess_i}.time = E_speechlocked_avg{sess_i}.time{1};
                    E_speechlocked_avg{sess_i}.trial = {squeeze(mean(cat(3,E_speechlocked{sess_i}.trial{:}),3,'omitnan'))};
                    % detrend trial average
                    %         E_cuelocked_detr{sess_i} =   E_cuelocked{sess_i};
                    %         E_cuelocked_detr{sess_i}.trial = cellfun(@(x) x - E_cuelocked_avg{sess_i}.trial{1,1}, E_cuelocked{sess_i}.trial,'uni',0);
                    E_speechlocked_detr{sess_i} =   E_speechlocked{sess_i};
                    E_speechlocked_detr{sess_i}.time = E_speechlocked_detr{sess_i}.time{1};
                    E_speechlocked_detr{sess_i}.trial = cellfun(@(x) x - E_speechlocked_avg{sess_i}.trial{1,1}, E_speechlocked{sess_i}.trial,'uni',0);

                    % calculate TF
                    disp(' TF decomposition in locked activity' )
                    cfg = [];
                    cfg.times = E_speechlocked_avg{sess_i}.time;
                    cfg.trials = 1;
                    cfg.pnts = numel(E_speechlocked_avg{sess_i}.time);
                    cfg.srate = E_speechlocked_avg{sess_i}.fsample;
                    cfg.trials = numel(E_speechlocked{sess_i}.trial);
                    E_speechlocked_avg{sess_i}.TF = wav_decomp(E_speechlocked_avg{sess_i}.trial{1,1},cfg,fmin,fmax,nfreq,'means',E_speechlocked_avg{sess_i}.time);




                    disp(' TF decomposition in (locked + non-locked) activity' )
                    cfg = [];
                    cfg.times = E_speechlocked{sess_i}.time{1,1};
                    cfg.trials = 1;
                    cfg.pnts = numel(E_speechlocked{sess_i}.time{1,1});
                    cfg.srate = E_speechlocked{sess_i}.fsample;
                    cfg.trials = numel(E_speechlocked{sess_i}.trial);
                    E_speechlocked{sess_i}.TF = wav_decomp(cat(3,E_speechlocked{sess_i}.trial{:}),cfg,fmin,fmax,nfreq,'means',E_speechlocked{sess_i}.time{1,1});
                    % calculate TF baseline norm
                    baselinetime = [(epochs.stim1_onset - .6) - epochs.syl1_onset  (epochs.stim1_onset - .1) - epochs.syl1_onset];
                    E_speechlocked{sess_i}.TFnorm = wav_decomp(E_speechlocked{sess_i}.trial{1,1},cfg,fmin,fmax,nfreq,'means',E_speechlocked{sess_i}.time{1,1},baselinetime);



                    disp(' TF decomposition in non-locked activity' )
                    cfg = [];
                    cfg.times = E_speechlocked_detr{sess_i}.time;
                    cfg.trials = 1;
                    cfg.pnts = numel(E_speechlocked_detr{sess_i}.time);
                    cfg.srate = E_speechlocked_detr{sess_i}.fsample;
                    cfg.trials = numel(E_speechlocked_detr{sess_i}.trial);
                    E_speechlocked_detr{sess_i}.TF = wav_decomp(cat(3,E_speechlocked_detr{sess_i}.trial{:}),cfg,fmin,fmax,nfreq,'means',E_speechlocked_detr{sess_i}.time);



            end

        end

    end

    % save current subject for Cue
    fprintf('Saving TF results for Subject %s %d - %d  in Event %s \n', SUBJECT, subj_i, n_SUBJECTS, events{1})
    PATH_TF = fullfile(PATHS.DataDir ,SUBJECT,'TF_speech.mat');
    save(PATH_TF, 'E_speechlocked','E_speechlocked_avg','E_speechlocked_detr','-v7.3')
    % save current subject for Speech
    fprintf('Saving TF results for Subject %s %d - %d  in Event %s \n', SUBJECT, subj_i, n_SUBJECTS, events{2})
    PATH_TF = fullfile(PATHS.DataDir ,SUBJECT,'TF_cue.mat');
    save(PATH_TF, 'E_cuelocked','E_cuelocked_avg','E_cuelocked_detr','-v7.3')
    
    % clear some variables
        
    clear E epochs E_cuelocked E_cuelocked_avg E_cuelocked_detr E_speechlocked E_speechlocked_avg E_speechlocked_detr


end

    %%





    %%
    figure('Renderer','painters','Position',[300 300 700 300])
    hold on
    plot(E_speechlocked_avg{1,1}.time, ((E_speechlocked_avg{1,1}.trial{1,1}(1,:))),'r')
    ylabel('ECoG amp. [a.u.]')
    axis("tight")
    box off
    xlabel(" Time w.r.t. speech onset")

    cfg = [];
    cfg.times = E_speechlocked_avg{1,1}.time;
    cfg.trials = 1;
    cfg.pnts = numel(E_speechlocked_avg{1,1}.time);
    cfg.srate = E_speechlocked_avg{1,1}.fsample;
    cfg.trials = numel(E_speechlocked{1,1}.trial);
    tf_locked = wav_decomp(cat(3,E_speechlocked{1,1}.trial{:}),cfg,2,140,60,'means',E_speechlocked_avg{1,1}.time);
    tf_locked_detr = wav_decomp(cat(3,E_speechlocked_detr{1,1}.trial{:}),cfg,fmin,fmax,nfreq,'means',E_speechlocked_avg{1,1}.time);


    figure('Renderer','painters','Position',[300 300 1200 300])
    tiledlayout(1,3)
    nexttile
    %     imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked_avg(10,:,:)).*repmat(logspace(log10(2),log10(140),60),numel(E_speechlocked_avg{1,1}.time),1)');
    imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked_avg(10,:,:)));

    colormap(linspecer)
    colorbar
    set(gca,'YDir','normal')
    hold on
    plot(E_speechlocked_avg{1,1}.time, 110 + zscore((E_speechlocked_avg{1,1}.trial{1,1}(1,:))),'color','w','LineWidth',1.5)
    clim([min([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)])  max([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)]) ])

    title('Speech-locked')
    nexttile
    % imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked(10,:,:)).*repmat(logspace(log10(2),log10(140),60),numel(E_speechlocked_avg{1,1}.time),1)');
    imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked(10,:,:)));

    colormap(linspecer)
    colorbar
    title('Speech-locked + Non-speech-locked')
    clim([min([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)])  max([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)]) ])

    nexttile
    imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked_detr(10,:,:)));
    %imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked_detr(10,:,:)).*repmat(logspace(log10(2),log10(140),60),numel(E_speechlocked_avg{1,1}.time),1)');
    colormap(linspecer)
    colorbar
    title(' Non-speech-locked')
    clim([min([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)])  max([tf_locked_detr(:); tf_locked(:); tf_locked_avg(:)]) ])




    figure('Renderer','painters','Position',[300 300 1200 300])

    imagesc(E_speechlocked_avg{1,1}.time, logspace(log10(2),log10(140),60),squeeze(tf_locked_detr(10,:,:) - tf_locked(10,:,:)).*repmat(logspace(log10(2),log10(140),60),numel(E_speechlocked_avg{1,1}.time),1)');
    colormap(linspecer)
    colorbar
    title(' Speech-locked (difference)')



    figure('Renderer','painters','Position',[300 300 700 300])
    plot(E_speechlocked_avg{1,1}.time, squeeze(cc(1,:,1:10)),'color','k','linewidth',.25)
    hold on
    plot(E_speechlocked_avg{1,1}.time, mean(squeeze(cc(1,:,1:end)),2,'omitnan')','color','r','linewidth',1)

    xlim([-2.5 2])
    for sess_i =  1 : nSessions
        cfg_ft = [];
        cfg_ft.channel = E_channel{ch_i};
        Edata = ft_selectdata(cfg_ft,E);
    end
    %
    % calculate Evoked activity


    %
    %     for ch_i = 1 : nEChannels
    %
    %
    %
    %     E_channel = E_channels{pairs_idx(pair_i,1)};
    %
    %
    %     for sess_i = 1 : nSessions
    %         E_epoch = E.epochs(E.epochs.session_id_coding == sess_i,:);
    %         x_ = [];
    %         for ee = 1 : height(E_epoch)
    %             start_ = E_epoch.starts(ee);
    %             stop_ = E_epoch.ends(ee);
    %             idx_ = E_data.time{sess_i} >= start_ & E.time{sess_i} <= stop_;
    %             x_ = [x_ ;E_data.trial{sess_i}(idx_)];
    %         end
    %         E_evoked.trial{sess_i}(idx_) = mean(x_,'omitnan');
    %         for ee = 1 : height(E_epoch)
    %             start_ = E_epoch.starts(ee);
    %             stop_ = E_epoch.ends(ee);
    %             idx_ = E_data.time{sess_i} >= start_ & E.time{sess_i} <= stop_;
    %             E_data.trial{sess_i}(idx_) = E_data.trial{sess_i}(idx_) - evoked;
    %         end
    %     end
    %     end



    %
    %         STS.sel_unit = sel_units;
    %         STS.n_units = numel(STS.label);
    %         STS.n_channels = numel(E.label);

    %     % save plv value
    %     disp("save plv results... ")
    %     mybeep
    %     save(fullfile(PATHS.saveMatFiles,SUBJECT,'PLV.mat'), 'PLV',"-v7.3")
    %
    % end
    %
    % tEnd = toc(tStart);
    % disp(" Script ends ...")
    % fprintf(" Time elapsed %.2f \n",tEnd)
    % mybeep % play sound of end script

