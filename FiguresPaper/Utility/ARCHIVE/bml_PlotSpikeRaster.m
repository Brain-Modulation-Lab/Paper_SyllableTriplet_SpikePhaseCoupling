function bml_PlotSpikeRaster(triplets_events,produced_syllables,stimulus_syllables,PATH_OUTPUT,PATH_FIGURES,SUBJECT,session_id,unit_id,sortnotes,event0_key,fs)

try
    n_event0 = numel(event0_key);
    figure("renderer","painters")
    sgtitle(strcat(SUBJECT,'_Session', num2str(session_id),'_Channel',sortnotes.channel{unit_id},'_Unit',num2str(sortnotes.unit(unit_id))),"interpreter","None")
    % lock units
    for event0_i = 1 : n_event0
        if contains(event0_key{event0_i},"cue")
            event0 = triplets_events.stim_onset;
            RT_rel = triplets_events.prod_onset - event0;
            [~,ord] = sort(RT_rel);
            c_event = "m";
            SpDur_rel = triplets_events.prod_offset - event0;
            respInterval = round(fs*[-nanmedian(RT_rel) - 0.5 nanmedian(SpDur_rel)+0.5])/fs;
        elseif contains(event0_key{event0_i},"stim")
            event0 = triplets_events.prod_onset;
            RT_rel = triplets_events.prod_onset - event0;
            StimStop_rel = triplets_events.stim_offset - event0;
            c_event = "g";
            
            SpDur_rel = triplets_events.prod_offset- event0;
            [~,ord] = sort(SpDur_rel);
            respInterval = round(fs*[nanmedian(StimStop_rel) - 0.5 nanmedian(SpDur_rel)+0.5])/fs;
        else
            warning("Please eneter a valid event key (cue or stim)")
        end
        
        
        PATH_IFRSTIM = fullfile(PATH_OUTPUT,strcat(SUBJECT,'_Session', num2str(session_id),'_Channel',sortnotes.channel{unit_id},'_Unit',num2str(sortnotes.unit(unit_id)),'_',event0_key{event0_i})); %s, Unit %d\n,'_IFR_Zstat'))
        %Nobs = round((respInterval(2)-respInterval(1))/(2*filtSD/fs));
        load(PATH_IFRSTIM);
        if contains(event0_key{event0_i},"cue")
            Stim = CueOnset;
            clear CueOnset
        elseif contains(event0_key{event0_i},"stim")
            Stim = SpeechOnset;
            clear SpeechOnset
        end
        tvect = Stim.respInterval(1):1/fs:Stim.respInterval(2);
        ntrials = numel(Stim.trial_samples);
        
        subplot(3,n_event0,event0_i)
        
        hold on
        cont = 1;
        for trial_i = ord'
            syllable = produced_syllables(produced_syllables.session_id == session_id & produced_syllables.trial_id == triplets_events.trial_id(trial_i),:);
            n_syllables = height(syllable);
            syllable_rel_start = syllable.starts - event0(triplets_events.trial_id == triplets_events.trial_id(trial_i));
            syllable_rel_end = syllable.ends - event0(triplets_events.trial_id == triplets_events.trial_id(trial_i));
            
            spkstimes = tvect(Stim.DD(trial_i,:)==1);
            nspks = numel(spkstimes);
            scatter(spkstimes,cont*ones(1,nspks),1.5,'k','filled')
            plot([syllable_rel_start syllable_rel_end]',cont*ones(2,n_syllables),'g')
            if ~ismissing(stimulus_syllables)
                stim_syllable = stimulus_syllables(stimulus_syllables.session_id == session_id & stimulus_syllables.trial_id == triplets_events.trial_id(trial_i),:);
                stim_syllable_rel_start = stim_syllable.starts - event0(triplets_events.trial_id == triplets_events.trial_id(trial_i));
                stim_syllable_rel_end = stim_syllable.ends - event0(triplets_events.trial_id == triplets_events.trial_id(trial_i));
                plot([stim_syllable_rel_start stim_syllable_rel_end]',cont*ones(2,n_syllables),'m')
            end
            cont = cont + 1;
        end
        ylabel(" Trial ID ")
        xlim(Stim.respInterval)
        meanIFRStim = mean(Stim.IFRdata);
        UpIFRStim = meanIFRStim + std(Stim.IFRdata)/sqrt(ntrials);
        DownIFRStim = meanIFRStim - std(Stim.IFRdata)/sqrt(ntrials);
        meanISITSStim = mean(Stim.ISITSdata);
        UpISITSStim = meanISITSStim + std(Stim.ISITSdata)/sqrt(ntrials);
        DownISITSStim = meanISITSStim - std(Stim.ISITSdata)/sqrt(ntrials);
        ax1 = gca;
        ax1.XAxis.Visible = 'off';
        line([0 0],[0 ntrials],"linewidth",2,"color",c_event)
        title(event0_key{event0_i})
        
        subplot(3,n_event0,event0_i + n_event0)
        plot(tvect,meanIFRStim,"linewidth",1.5,"color","k")
        hold on
        fill([tvect fliplr(tvect)],  [UpIFRStim fliplr(DownIFRStim)], 'k',"facealpha",.15);
        line([0 0],[0 ntrials],"linewidth",2,"color",c_event)
        set(gca,"box","off")
        xlabel(" Time [s] ")
        ylabel(" Firing Rate [Hz]")
        grid on
        grid minor
        signclusters_idx = find([Stim.IFRmod.Zflag] == 1 & [Stim.IFRmod.Zlength_enough] == 1);
        nsignclusters = numel(signclusters_idx);
        for clust_i = 1 : nsignclusters
            if sign(Stim.IFRmod(signclusters_idx(clust_i)).Zclust) > 0
                plot(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins),max(meanIFRStim + 2.5*std(Stim.IFRdata)/sqrt(ntrials))*ones(size(Stim.IFRmod(signclusters_idx(clust_i)).tbins)),"r","linewidth",4)
                plot(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins),meanIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins),"linewidth",1.5,"color","r")
                fill([tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins) fliplr(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins))],  [UpIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins) fliplr(DownIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins))], 'r',"facealpha",.25);
            else
                plot(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins),max(meanIFRStim + 2.5*std(Stim.IFRdata)/sqrt(ntrials))*ones(size(Stim.IFRmod(signclusters_idx(clust_i)).tbins)),"b","linewidth",4)
                plot(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins),meanIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins),"linewidth",1.5,"color","b")
                fill([tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins) fliplr(tvect(Stim.IFRmod(signclusters_idx(clust_i)).tbins))],  [UpIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins) fliplr(DownIFRStim(Stim.IFRmod(signclusters_idx(clust_i)).tbins))], 'b',"facealpha",.25);
            end
        end
        ylim([min(meanIFRStim - 2*std(Stim.IFRdata)/sqrt(ntrials)) max(meanIFRStim + 3*std(Stim.IFRdata)/sqrt(ntrials))]);
        IFRbase_CI = [median(meanIFRStim) - 1.5*iqr(meanIFRStim) median(meanIFRStim) + 1.5*iqr(meanIFRStim)];;
        plot(Stim.respInterval, [IFRbase_CI(1) IFRbase_CI(1)],"--k","linewidth",1.5)
        plot(Stim.respInterval, [IFRbase_CI(2) IFRbase_CI(2)],"--k","linewidth",1.5)
        xlim(Stim.respInterval)
        
        subplot(3,n_event0,event0_i + 2*n_event0)
        plot(tvect,meanISITSStim,"linewidth",1.5,"color","k")
        hold on
        fill([tvect fliplr(tvect)],  [UpISITSStim fliplr(DownISITSStim)], 'k',"facealpha",.15);
        line([0 0],[0 ntrials],"linewidth",2,"color",c_event)
        set(gca,"box","off")
        xlabel(" Time [s] ")
        ylabel(" ISI [s]")
        grid on
        grid minor
        signclusters_idx = find([Stim.ISITSmod.Zflag] == 1 & [Stim.ISITSmod.Zlength_enough] == 1);
        nsignclusters = numel(signclusters_idx);
        for clust_i = 1 : nsignclusters
            if sign(Stim.ISITSmod(signclusters_idx(clust_i)).Zclust) > 0
                plot(tvect(Stim.ISITSmod(signclusters_idx(clust_i)).tbins),max(meanISITSStim + 2.5*std(Stim.ISITSdata)/sqrt(ntrials))*ones(size(Stim.ISITSmod(signclusters_idx(clust_i)).tbins)),"b","linewidth",4)
                plot(tvect(Stim.ISITSmod(signclusters_idx(clust_i)).tbins),meanISITSStim(Stim.ISITSmod(signclusters_idx(clust_i)).tbins),"linewidth",1.5,"color","b")
                fill([tvect(Stim.ISITSmod(signclusters_idx(clust_i)).tbins) fliplr(tvect(Stim.ISITSmod(signclusters_idx(clust_i)).tbins))],  [UpISITSStim(Stim.ISITSmod(signclusters_idx(clust_i)).tbins) fliplr(DownISITSStim(Stim.ISITSmod(signclusters_idx(clust_i)).tbins))], 'b',"facealpha",.25);
            end
        end
        ylim([min(meanISITSStim - 2*std(Stim.ISITSdata)/sqrt(ntrials)) max(meanISITSStim + 3*std(Stim.ISITSdata)/sqrt(ntrials))]);
        
        ISITSbase_CI = [median(meanISITSStim) - 1.5*iqr(meanISITSStim) median(meanISITSStim) + 1.5*iqr(meanISITSStim)];
        plot(Stim.respInterval, [ISITSbase_CI(1) ISITSbase_CI(1)],"--k","linewidth",1.5)
        plot(Stim.respInterval, [ISITSbase_CI(2) ISITSbase_CI(2)],"--k","linewidth",1.5)
        xlim(Stim.respInterval)
        
        clear Stim
    end
    
    figname = fullfile(PATH_FIGURES,strcat(SUBJECT,'_Session', num2str(session_id),'_Channel',sortnotes.channel{unit_id},'_Unit',num2str(sortnotes.unit(unit_id)),'_locked',event0_key{event0_i})); %s, Unit %d\n,'_IFR_Zstat'))
    disp("displaying and saving figure...")
    saveas(gcf,figname,"png")
    saveas(gcf,figname,"pdf")
    saveas(gcf,figname,"eps")
    close gcf
catch 
    warning(" No enough samples")
end
    