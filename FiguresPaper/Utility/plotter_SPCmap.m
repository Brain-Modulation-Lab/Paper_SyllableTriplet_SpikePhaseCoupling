function fh = plotter_SPCmap(SPC, cfg)

% unfold cfg
Position  = cfg.Position;
EvtTimes = cfg.EvtTimes;
EvtTypes = cfg.EvtTypes;
if isfield(cfg,'CLim')
    CLim = cfg.CLim;
else
    CLim = [];
end

if isfield(cfg,'CMap')
    CMap = cfg.CMap;
end




fh{1} = figure('renderer','painters','position',Position);
imagesc(SPC.timeQ,SPC.freqQ,SPC.map)
hold on
for evt = [2 3 4 5]'
    plot(EvtTimes(evt)' *[1,1], [min(SPC.freqQ) max(SPC.freqQ)], '--', 'Color', 'w', 'LineWidth', 2)
end
set(gca, 'XTick', EvtTimes([2 3 4 5]))
xLabeling = EvtTypes(2:end-1);
set(gca, 'XTickLabel', xLabeling)
set(gca,'YDir','Normal')
xtickangle(30)
yticks(round(SPC.freqQ(1:6:end),2))
set(gca,'YMinorTick','off')
cb = colorbar;
ylabel(cb, 'PPC [z-score]')
ylabel('Frequency [Hz] ')
set(gca,'YScale','log')
colormap(CMap)
if ~isempty(CLim)
    clim(CLim)
end
axis tight
box off


if isfield(SPC,'Stat')
    fh{2} = figure('renderer','painters','position',Position);
    imagesc(SPC.timeQ,SPC.freqQ,SPC.Stat.TscoreMap)
    hold on
    for evt = [2 3 4 5]'
        plot(EvtTimes(evt)' *[1,1], [min(SPC.freqQ) max(SPC.freqQ)], '--', 'Color', 'w', 'LineWidth', 2)
    end
    set(gca, 'XTick', EvtTimes([2 3 4 5]))
    xLabeling = EvtTypes(2:end-1);
    set(gca, 'XTickLabel', xLabeling)
    set(gca,'YDir','Normal')
    xtickangle(30)
    yticks(round(SPC.freqQ(1:6:end),2))
    set(gca,'YMinorTick','off')
    cb = colorbar;
    ylabel(cb, 'PPC [z-score]')
    ylabel('Frequency [Hz] ')
    set(gca,'YScale','log')
    colormap(CMap)
%     if ~isempty(CLim)
%         clim(CLim)
%     end
    axis tight
    box off
    hold on
    [~,idx] = sort(SPC.Stat.Zstat,'descend');
    
    if isfield(cfg,'NClustplot')
        if isnumeric(cfg.NClustplot)
            NClustplot = cfg.NClustplot;
        else
            NClustplot = max(idx);
        end
    else
        NClustplot = 3;
    end

    for idx_i = 1 : NClustplot
    
            contour(SPC.timeQ,SPC.freqQ,SPC.Stat.LabelMatrix.*(SPC.Stat.LabelMatrix == idx(idx_i)),1, 'linecolor','k', 'linewidth', 2)
        
        %caxis([-10 10])
    end
end

