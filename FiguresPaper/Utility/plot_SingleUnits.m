function fh = plot_SingleUnits(cfg,DB)

% get units pairs
[Units, idxUnitPairs] = getUniquePairs(DB.Pairs.Location);
nUnits = size(Units,1);
nUniquePairs = accumarray(idxUnitPairs,1);

% get unit clusters
cfg.VarFields = {'Sign_Taskid','ClusterCentroidBand','ClusterSChannel','ClusterSubjects','ClusterS_FRMod','ClusterOnOff','ClusterCentroid'};
[Sign_Task,ClusterCentroidBand,Cluster_SUnit, Cluster_Subject, ClusterS_FRMod,ClusterOnOff,ClusterCentroid] = get_SPCClustersProperties(cfg,DB.Clusters);
%Units_FRmod = arrayfun(@(x) , unique(idxUnitPairs));


%[UnitsClusters, IdxUnitClusters] = getUniquePairs(DB.Clusters);
SPCbeta = nan(4,nUnits);
UnitFRMod = nan(1,nUnits);
SPCbandsunits_Density = table();
for unit_i = 1 : nUnits
    idxU = all(contains([Cluster_SUnit string(Cluster_Subject)],Units(unit_i,:)),2);
    sum(idxU)
    UnitClusterCentroidBand = ClusterCentroidBand(idxU);
    UnitClusterCentroidFrequency = ClusterCentroid(idxU,2);
    UnitClusterOnOff = ClusterOnOff(idxU,:)';
    UnitSign_Task = Sign_Task(:,idxU);
    %tmp = ClusterS_FRMod(idxU == 1);
    UnitFRMod(unit_i) = mean(DB.Pairs.Location.S_typeFRmod(idxUnitPairs == unit_i));
    SPCbeta(4,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(5,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(3,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(2,unit_i) = sum(UnitClusterCentroidBand == 3 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(1,unit_i) = sum(UnitClusterCentroidBand == 2 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);

    SPCbandsunits_Density.Unit(unit_i) = Units(unit_i,2);
    SPCbandsunits_Density.Subject(unit_i)  = Units(unit_i,1);
    SPCbandsunits_Density.nPairs(unit_i)  = nUniquePairs(unit_i);
    SPCbandsunits_Density.FRmod(unit_i) = UnitFRMod(unit_i);
    SPCbandsunits_Density.Density(unit_i)  = numel(UnitClusterCentroidBand)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.theta(unit_i)  = sum(UnitClusterCentroidBand == 2)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.alpha(unit_i)  = sum(UnitClusterCentroidBand == 3)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.beta(unit_i)  = sum(UnitClusterCentroidBand == 4)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.gammaL(unit_i)  = sum(UnitClusterCentroidBand == 5)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.gammaH(unit_i)  = sum(UnitClusterCentroidBand == 6)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.thetaalpha_ITI(unit_i)  = sum(UnitClusterCentroidBand < 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.thetaalpha_speech(unit_i)  = sum(UnitClusterCentroidBand < 4 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
        %SPCbandsunits_Density.beta_ITI(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1  & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;

   SPCbandsunits_Density.beta_ITI(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
  %SPCbandsunits_Density.beta_speech(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= -0.2 & UnitClusterOnOff(2,:) <= 0.5)/nUniquePairs(unit_i)*100;
  SPCbandsunits_Density.beta_speech(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitClusterOnOff(1,:) >= -0.2 & UnitClusterOnOff(2,:) <= 0.5)/nUniquePairs(unit_i)*100;
   
  % SPCbandsunits_Density.beta_speech(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(2,:) <= 0.5)/nUniquePairs(unit_i)*100;
      SPCbandsunits_Density.beta_postspeech(unit_i)  = sum(UnitClusterCentroidBand == 4  & UnitClusterOnOff(1,:) >= 0.51 & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
     % SPCbandsunits_Density.beta_postspeech(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(5,:) == 1 & UnitClusterOnOff(1,:) >= 0.51 & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    
    SPCbandsunits_Density.meanfreq(unit_i)  = mean(UnitClusterCentroidFrequency);
    SPCbandsunits_Density.medianfreq(unit_i)  = median(UnitClusterCentroidFrequency);

end
    
    
    
    
    %   SPCbeta(2,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(5,:) == 1)/nUniquePairs(unit_i);
%    SPCbeta(1,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1)/nUniquePairs(unit_i);





SPCbeta_unittypes = arrayfun(@(x) SPCbeta(:,UnitFRMod == x),[-1 0 1 2],'uni',false);
% (no, no),  (yes,no),(no, yes), (yes, yes)
type_cells = nan(4,4);
UnitTypes = {'D','N','I','M'};
for ui = 1 : 4
    xxx = SPCbeta_unittypes{ui};
    type_cells(ui,1) = sum(xxx(3,:) == 0 & xxx(4,:) == 0);
    type_cells(ui,2) = sum(xxx(3,:) > 0 & xxx(4,:) == 0);
    type_cells(ui,3) = sum(xxx(3,:) == 0 & xxx(4,:) > 0);
    type_cells(ui,4) = sum(xxx(3,:) > 0 & xxx(4,:) > 0);
    
    % save csv
    TableToSave = table();
    TableToSave.ITI = xxx(1,:)';
    TableToSave.PostSpeech = xxx(2,:)';

    writetable(TableToSave,['SPCbeta_unit-',UnitTypes{ui},'.csv'])
    
end

% start figures

fh{1} = figure('Position',[200 200  800 1200]);
nexttile
imagesc(SPCbeta')
colormap(linspecer)
colorbar

% for ui = 1 : 4
%     nexttile
%     imagesc( SPCbeta_unittypes{ui}')
%     %colormap(linspecer)
%     colormap(linspecer)
%     colorbar
%     clim([0 0.75])
%     
% end
nexttile
for ui = 1 : 4
    nexttile
    plot( [1;3] + 0.1*randn(2, size(SPCbeta_unittypes{ui},2)), (100*(SPCbeta_unittypes{ui}(3:4,:))+eps),...
        'linestyle','-','markerfacecolor',[.6 .6 .6],'markeredgecolor',[0 0 0], ...
        'color',[0 0 0],'markersize',8,'marker','o')
    %colormap(linspecer)
    ylim([0.01 10^1.5])
%     xlim([0.5 3.5])
    box off
    xticks([1 3])
    xticklabels({'ITI','Post-Speech'})
    ylabel(' SPC beta')
    set(gca,'yscale','log')
end



fh{2} = figure('Position',[200 200  800 800]);
for ui = 1 : 4
    nexttile
    scatter((100*SPCbeta_unittypes{ui}(3,:) + eps), ...
        (100*SPCbeta_unittypes{ui}(4,:) + eps), ...
        'markerfacecolor',[.8 .8 .8],'markeredgecolor',[.6 .6 .6], ...
        'color',[0 0 0],'SizeData',55,'marker','o')
    %colormap(linspecer)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
%     ylim([0 52])
%     xlim([0 52])
%     xticks(0:10:50)
%     yticks(0:10:50)
% ylim([10^-15 3*10^5])
% xlim([10^-15 3*10^5])

    box off
%     xticks([1 3])
%     xticklabels({'ITI','Post-Speech'})
    ylabel(' SPC beta')
    %set(gca,'yscale','log')
end


% figure
% for ui = 1 : 4
%     nexttile
% histogram(-diff(SPCbeta_unittypes{ui}),100,'normalization','probability')
% %colormap(linspecer) 
% % ylim([0 0.5])
% xlim([-0.2 0.5])
% end
% 
% 
% figure
% scatter(SPCbeta(1,:), SPCbeta(2,:))
% xlim([0 0.6])
% ylim([0 0.6])
% figure
% for ui = 1 : 4
%     nexttile
% scatter(SPCbeta_unittypes{ui}(1,:), SPCbeta_unittypes{ui}(2,:),'filled')
% colormap(linspecer)    
% end
% 
%markers = ['s','d','o','t']
% fh{2} = figure('Position',[200 200  800 800]);
% for ui = 1 : 4
%     nexttile
%     scatter(SPCbeta_unittypes{ui}(1,:), SPCbeta_unittypes{ui}(2,:),25,'filled')
%     colormap(linspecer)
%     %[r,p] = corr(SPCbeta_unittypes{ui}(1,:)', SPCbeta_unittypes{ui}(2,:)')
%     hold on
%     plot([0 0.6], [0 0.6], '--r')
%     xlim([0 0.6])
%     ylim([0 0.6])
% end



% SU_cluster = table();
% cnt_tbl = 1;
% 
% % generate surrogates of random spectral specificity
% % btsp = 1000;
% ClusterSubjects_unittype = arrayfun(@(x) ClusterSubjects(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
% ClusterSChannel_unittype = arrayfun(@(x) ClusterSChannel(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
% SignTask_id_unittype = arrayfun(@(x) SignTask_id(:,ClusterS_FRMod == x), unit_types,'UniformOutput',false);
% cluster_prob_unittypes = arrayfun(@(x) cluster_prob(ClusterS_FRMod == x, : ), unit_types,'UniformOutput',false);
% ClusterCentroidBand_unittypes = arrayfun(@(x) ClusterCentroidBand(ClusterS_FRMod == x), unit_types,'UniformOutput',false);
% 
% 
% 
% ClusterSingleUnitsSequence = cell(1,numel(unit_types));
% check = [];
% for unit_i = 1 : numel(unit_types)
% 
%     hold on
%     [uniqueUnits,idauniqueUnits,idxuniqueUnits] = unique([ClusterSubjects_unittype{unit_i} ClusterSChannel_unittype{unit_i}],'rows');
%     nuniqueUnits_clus = size(uniqueUnits,1);
%     uniqueSubjects = unique(uniqueUnits(:,1));
%     nuniqueSubjects = numel(uniqueSubjects);
% 
% 
%     nClusters_uniqueUnits = accumarray(idxuniqueUnits,1);
%     modefrequencybandsUnits  = arrayfun(@(x) mode(ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == x)),unique(idxuniqueUnits));
%     [modefrequencybandsUnits_sorted, idx_modefrequencybandsUnits_sorted] = sort(modefrequencybandsUnits);
%     idxuniqueUnits_sorted = idxuniqueUnits(idx_modefrequencybandsUnits_sorted);
% 
%     %uniqueUnitstrans = [1 ;find(abs(diff(idxuniqueUnits)) > 0) + 1];
%     cnt = nClusters_pertype(unit_i);
%     %yax = linspace(nClusters_pertype(unit_i), 0,nClusters_pertype(unit_i));
%     for clust_i = idx_modefrequencybandsUnits_sorted'
%         SU_cluster.SubjID(cnt_tbl) =  uniqueUnits(clust_i,1);
%         SU_cluster.UnitID(cnt_tbl) =  uniqueUnits(clust_i,2);
%         SU_cluster.UnitTypeFR(cnt_tbl) =  unit_types(unit_i);
% 
%         tmp_clusters = cluster_prob_unittypes{unit_i}((idxuniqueUnits == clust_i),:);
%         ClusterSingleUnitsSequence{unit_i}{end+1} = tmp_clusters;
%         tmp_signtask = SignTask_id_unittype{unit_i}(:,idxuniqueUnits == clust_i);
%         tmp_bands = ClusterCentroidBand_unittypes{unit_i}(idxuniqueUnits == clust_i);
%         SU_cluster.nClusters(cnt_tbl) =  numel(tmp_bands);
%         check = [check; [SU_cluster.nClusters(cnt_tbl) sum(contains(ClusterSubjects,SU_cluster.SubjID(cnt_tbl)) & strcmpi(ClusterSChannel,SU_cluster.UnitID(cnt_tbl))) ]];
%         clusters_probands_byunit = arrayfun(@(x) mean(tmp_bands == x), 1:6);
% 
%         for fi = 2:6
%             SU_cluster.(char(bands.name(fi)))(cnt_tbl) = clusters_probands_byunit(fi);
%         end
% %         Entropyband = sum(-clusters_probands_byunit.*log2(clusters_probands_byunit + eps));
% %         if Entropyband < 0
% %             Entropyband = 0;
% %         elseif Entropyband > 1
% %             Entropyband = 1;
% %         end
% % 
% %         SU_cluster.BandSpecifity(cnt_tbl) = 1 - Entropyband/log2(6);
% %         BandSpecifity_btsp = nan(1,btsp);
% %         btsp_values = randi(5,btsp,numel(tmp_bands));
% %         btsp_values_p = nan(btsp,6);
% %         for btsp_i = 1 : btsp
% %             btsp_values_p(btsp_i,:) = arrayfun(@(x) mean(btsp_values(btsp_i,:) == x,2),1:6);
% %             BandSpecifity_btsp(btsp_i) = 1 - sum(-btsp_values_p(btsp_i,:).*log2(btsp_values_p(btsp_i,:) + eps))/log2(6);
% %         end
% %         SU_cluster.BandSpecifity_prct95(cnt_tbl) = prctile(BandSpecifity_btsp,95);
% %         SU_cluster.BandSpecifity_corr(cnt_tbl) = 1 - (Entropyband - mean(BandSpecifity_btsp,'omitnan')) /(log2(6) - mean(BandSpecifity_btsp,'omitnan')) ;
% %         SU_cluster.BandSpecifity_prct95_corr(cnt_tbl) = 1 - (SU_cluster.BandSpecifity_prct95(cnt_tbl)- mean(BandSpecifity_btsp,'omitnan')) /(log2(6) - mean(BandSpecifity_btsp,'omitnan')) ;
% 
%         SU_cluster.nPairs(cnt_tbl) =  sum(contains(PairSubjects_unittype{unit_i}, uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),1)) & contains(PairSChannel_unittype{unit_i},uniqueUnits(idx_modefrequencybandsUnits_sorted(clust_i),2)));
%         SU_cluster.ClusterDensity(cnt_tbl) =  SU_cluster.nClusters(cnt_tbl)/SU_cluster.nPairs(cnt_tbl);
% 
%         for fi = 2:6
%             for evt_i = 1 : 5
%                 tmp = sum(tmp_signtask(:,tmp_bands == fi),2)/SU_cluster.nPairs(cnt_tbl);
%                 SU_cluster.([char(bands.name(fi)),'_', char(Task_labels{evt_i})])(cnt_tbl) = tmp(evt_i);
%             end
%         end
% 
%         % need to divide by task segment
% 
%         cnt_tbl = cnt_tbl + 1;
% 
%     end
% 
% end
% 
% % SU_cluster.BandSpecifity_corr(SU_cluster.BandSpecifity_corr < 0) = 0;
% % SU_cluster.BandSpecifity_corr(SU_cluster.BandSpecifity_corr > 1) = 1;
% 
% % now plot the table
% matrix_clusterbands = cell(1,numel(unit_types));
% SU_cluster_tmp = arrayfun(@(x) SU_cluster(SU_cluster.UnitTypeFR == x,:), unit_types,'uni',false);
