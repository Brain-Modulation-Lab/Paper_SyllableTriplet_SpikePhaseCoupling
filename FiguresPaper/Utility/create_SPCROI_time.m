function [fh, ROI_AtlasLabels] = create_SPCROI_time(cfg, DB)

atlas = cfg.atlas;
min_subj = cfg.min_subj;
min_pairs = cfg.min_pairs;
flag_sort = cfg.flag_sort;
T = cfg.T;
nperms = cfg.nperms;
EvtTimes = cfg.EvtTimes;
EvtTypes = cfg.cfg.plot.EvtTypes;
bands = cfg.bands;

Pairs = DB.Pairs;
Clusters = DB.Clusters;

% get ROIs
ECoG_atlas = Pairs.Location.(atlas);
ECoG_subject = cellstr(Pairs.Location.subj_id);
[ROI_AtlasLabels,E_ROIAtlasLabel_pairs_n]  = getROIs(ECoG_atlas, ECoG_subject, min_subj,min_pairs,flag_sort);

% get variables
cfg.VarFields = {'clust_prob','ClusterEAtlasLabel','ClusterCentroidBand'};
[cluster_prob,ClusterEAtlasLabel, ClusterCentroidBand] = get_SPCClustersProperties(cfg, Clusters);


% analyze it
clust_prob_ROIs_bands = nan(6,numel(ROI_AtlasLabels),numel(T));

for fi = 2:6
    for roi_i = 1:numel(ROI_AtlasLabels)
        clust_prob_ROIs_bands(fi,roi_i,:) =  100*sum(cluster_prob(contains(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i})' & ClusterCentroidBand == fi,:),1)/sum(contains(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i})) ;
    end
end

% fh = figure('renderer','painters','position',[300 10 600 1600]);
% tiledlayout(5,2)
%clust_prob_ROIs_list_all = cell(1,6);
clust_prob_ROIs_sign_all = cell(1,6);
clust_prob_ROIs_perm_all = cell(1,6);
for fi = 2:6

    % perform clustering
    clust_prob_ROIs_list = cellfun(@(x) cluster_prob(contains(ClusterEAtlasLabel,x)' & ClusterCentroidBand == fi,:),ROI_AtlasLabels,'UniformOutput',false);


    clust_prob_ROIs_perm = nan(numel(ROI_AtlasLabels),numel(T),nperms);
    for roi_i  = 1 : numel(ROI_AtlasLabels)
        %clust_prob_ROI_perms = nan(nperms,numel(T));
        for perm_i = 1 : nperms
            tmp = [];
            cuts = randi(numel(T),1, size(clust_prob_ROIs_list{roi_i},1));
            for cut_i = 1 : numel(cuts)
                tmp = [tmp; clust_prob_ROIs_list{roi_i}(cut_i,[cuts(cut_i) : numel(T) 1 : (cuts(cut_i)-1)])];
            end
            clust_prob_ROIs_perm(roi_i,:,perm_i) = 100*sum(tmp)/sum(contains(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i}));
        end
    end
    %
    clust_prob_ROIs_lthr = prctile(clust_prob_ROIs_perm,5,3);
    clust_prob_ROIs_uthr = prctile(clust_prob_ROIs_perm,95,3);

    clust_prob_ROIs_sign = zeros(size(clust_prob_ROIs_bands,2),size(clust_prob_ROIs_bands,3));
    clust_prob_ROIs_sign(bsxfun(@lt,squeeze(clust_prob_ROIs_bands(fi,:,:)),clust_prob_ROIs_lthr)) = -1;
    clust_prob_ROIs_sign(bsxfun(@gt,squeeze(clust_prob_ROIs_bands(fi,:,:)),clust_prob_ROIs_uthr)) = 1;
    clust_prob_ROIs_sign(isnan(squeeze(clust_prob_ROIs_bands(fi,:,:)))) = nan;
    %     % perform plot
    %     nexttile(1+2*(fi-2),[1 2])
    %     imagesc(T,1:numel(ROI_AtlasLabels),squeeze(clust_prob_ROIs_bands(fi,:,:)))
    %     hold on
    %     % hold on
    %     for ii = 1: numel(ROI_AtlasLabels)
    %         scatter(T(clust_prob_ROIs_sign(ii,:) == 1 ),ii*ones(1,sum(clust_prob_ROIs_sign(ii,:) == 1)),5,'*r')
    %         scatter(T(clust_prob_ROIs_sign(ii,:) == -1 ),ii*ones(1,sum(clust_prob_ROIs_sign(ii,:) == -1)),5,'*b')
    %         scatter(T(isnan(clust_prob_ROIs_sign(ii,:))),ii*ones(1,sum(isnan(clust_prob_ROIs_sign(ii,:)))),5,'*k')
    %
    %         %scatter(T(store{ii}.LabelMatrix(ii,:) ~= 0 ),ii*ones(1,sum(store.LabelMatrix(ii,:) ~= 0)),5,'*k')
    %     end
    %     colormap(linspecer)
    %     yticks(1:numel(ROI_AtlasLabels));
    %     yticklabels(ROI_AtlasLabels)
    %     set(gca,'TickLabelInterpreter','none')
    %     cb = colorbar;
    %     ylabel(cb,'Overlap Cluster [%]')
    %     box off
    %     for evt = [2 3 4 5]'
    %         plot(EvtTimes(evt)' *[1,1], [1 numel(ROI_AtlasLabels)], '--', 'Color', 'w', 'LineWidth', 1)
    %     end
    %     set(gca, 'XTick', EvtTimes([2 3 4 5]))
    %     xLabeling = EvtTypes(2:end-1);
    %     set(gca, 'XTickLabel', xLabeling)
    %     xtickangle(30)
    %     caxis([0 17])
    %     title(bands.symbol(fi))
    %     %
    %clust_prob_ROIs_list_all{1,fi} = clust_prob_ROIs_list;
    clust_prob_ROIs_sign_all{1,fi} = clust_prob_ROIs_sign;
    clust_prob_ROIs_perm_all{1,fi} = clust_prob_ROIs_perm;
end

%

fh = cell(1, numel(ROI_AtlasLabels));
for roi_i = 1:numel(ROI_AtlasLabels)
    fh{roi_i} = figure('renderer','painters','position',[300 10 600 1600],'Visible','off');
    tiledlayout(5,1)

    for fi = 2:6
        nexttile
        hold on
        toplot  = squeeze(clust_prob_ROIs_bands(fi,roi_i,:));
        plot(T,toplot,'linewidth',1,'color',hex2rgb(bands.color(fi)));

        plotshaded(T,[squeeze(prctile(clust_prob_ROIs_perm_all{fi}(roi_i,:,:),5,3)); ...
            squeeze(prctile(clust_prob_ROIs_perm_all{fi}(roi_i,:,:),50,3)); ...
            squeeze(prctile(clust_prob_ROIs_perm_all{fi}(roi_i,:,:),95,3))], [0 0 0], .2)
        switch fi
            case 2
                title(ROI_AtlasLabels{roi_i})
                ylim([0 3])
            case 3
                ylim([0 13])
            case 4
                ylim([0 18])
            case 5
                ylim([0 9])
            case 6
                ylim([0 9])
        end
        box off
    end
    
    %legend(bands.symbol(2:6))
end
% figname = fullfile(PATHS.saveFigures, 'group-level',['timepercsign_perROI_byband',atlas]);
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')
% saveas(gcf,figname,'fig')