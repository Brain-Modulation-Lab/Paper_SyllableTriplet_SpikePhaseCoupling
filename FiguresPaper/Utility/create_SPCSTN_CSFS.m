function fh = create_SPCSTN_CSFS(cfg, DB);
EvtTypes = cfg.cfg.plot.EvtTypes;
EvtTimes = cfg.EvtTimes;
targets = cfg.targets;
T = cfg.T;
nperms = cfg.nperms;
atlas = cfg.atlas;
min_subj = cfg.min_subj;
min_pairs = cfg.min_pairs;
flag_sort = cfg.flag_sort;
bands = cfg.bands;
% load clusters
SPCDensity = readtable(cfg.map);
% load stn
load(cfg.ref,'atlases'); % manually load definition of DISTAL atlas.

keep_labels = {'STN.','STN_motor','STN_associative','STN_limbic'};
idx_labels = cellfun(@(x) find(contains(atlases.names,x)), keep_labels);

STN=arrayfun(@(x) reducepatch(atlases.roi{x,2}.fv,.5),idx_labels,'UniformOutput',false); % extract the left side


ECoG_atlas = DB.Pairs.Location.(atlas);
ECoG_subject = cellstr(DB.Pairs.Location.subj_id);
[ROI_AtlasLabels,ROI_nPairs]  = getROIs(ECoG_atlas, ECoG_subject, min_subj,min_pairs,flag_sort);
% remove some areas like MTG, pars. O trian and STG plan
ROI_AtlasLabels(end-2:end) = [];
ROI_nPairs(end-2:end) = [];   


% find coordinates (Y for CS and Z for FS)

cfg.VarFields = {'ClusterPhase','ClusterCentroid','ClusterCentroidBand','ClusterS_MNI','ClusterEAtlasLabel','ClusterDur','clust_prob'};
[ClusterPhase,ClusterCentroid,ClusterCentroidBand, ClusterS_MNI,ClusterEAtlasLabel,ClusterDur,cluster_prob] = get_SPCClustersProperties(cfg,DB.Clusters);

% average and peaks of SPC for each frequency band
ClusterS_MNI_Band_avg = arrayfun(@(x) mean(ClusterS_MNI(ClusterCentroidBand == x,:),'omitnan'),2:6,'uni',false);
ClusterS_MNI_Band_peak = nan(5,3);

for fi = 1:5
    %tmp = SPCDensity(SPCDensity.Color == fi,:);
    %[~,idx_] = max(tmp.Size,[],'omitnan');
    %ClusterS_MNI_Band_peak(fi,:) = tmp{idx_,{'X','Y','Z'}};
    tmp = SPCDensity{:,fi+4};
    [~,idx_] = max(tmp,[],'omitnan');
    ClusterS_MNI_Band_peak(fi,:) = SPCDensity{idx_,{'X_MNI','Y_MNI','Z_MNI'}};
end






% % check focality
% 
% figure('position',[200 200 1800 300])
% tiledlayout(1,5)
% 
% for fi = 1:5
%     coords = ClusterS_MNI(ClusterCentroidBand == fi+1,:);
%     true_dist_clust = mean(vecnorm(coords - mean(coords,'omitnan'),2,2),'omitnan'); % mm
%     perm_dist_clust = nan(1,nperms);
%     for perm_i = 1 : nperms
%         perm_loc = ClusterS_MNI(randi(height(ClusterS_MNI),1,size(coords,1)),:);
%         perm_dist_clust(perm_i) = mean(vecnorm(perm_loc - mean(perm_loc,'omitnan'),2,2),'omitnan'); % mm
%     end
%     nexttile
%     histogram(perm_dist_clust,round(sqrt(nperms)),'FaceColor',[.8 ,.8 ,.8])
%     hold on
%     box off
%     xlabel('Distance to centre [mm]')
%     ylabel(' # observations ')
%     xlim([0.5 3.5])
%     if true_dist_clust > min(perm_dist_clust)
%         p_val = (sum((perm_dist_clust - mean(perm_dist_clust)) <= (true_dist_clust - mean(perm_dist_clust)))+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
%     else
%         p_val = 1/nperms;
%     end
%     xline(true_dist_clust,'label',{sprintf('p = %1.3f',p_val)},'Color','r','linewidth',1.4,'LabelOrientation','horizontal','FontSize',15)
% end

%% coordinates plot

CoordBands = cell(1,6);
for fi = 2:6
    CoordBands{fi} = ClusterS_MNI(ClusterCentroidBand == fi,:) - mean(STN{1}.vertices);
end

pairs_test = nchoosek(2:6,2);
p_test = cell(1,3);
diff_test = cell(1,3);
for xi = 1 : 3
    p_test{xi} = nan(size(pairs_test,1),1);
    diff_test{xi} =  nan(size(pairs_test,1),1);
    for pi = 1 :  size(pairs_test,1)
        [p_test{xi}(pi), diff_test{xi}(pi)]  = permutation_diffTest_indsamples(CoordBands{pairs_test(pi,1)}(:,xi), CoordBands{pairs_test(pi,2)}(:,xi), 500);
    end
end

figure('position',[200 200 900 600])
tiledlayout(5,3)
for fi = 2: 6
    nexttile
    raincloud_plot(CoordBands{fi}(:,1), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,1)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,1))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

        nexttile
    raincloud_plot(CoordBands{fi}(:,2), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,2)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

        nexttile
    raincloud_plot(CoordBands{fi}(:,3), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,3)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,3))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

end


%%

STN_vx = STN{1}.vertices;
[coeff,score,latent,~,explained] = pca(zscore(STN_vx));

% how to apply transform (((STN_PCtbl{1,1:3} - mean(STN_vx))./std(STN_vx))*coeff).*std(STN_vx)

%scale_fac_dots = 5;
scale_fac_pc = 5;

STN_PCtbl = table(STN_vx(:,1),STN_vx(:,2),STN_vx(:,3),score(:,1)*std(STN_vx(:,1)),score(:,2)*std(STN_vx(:,2)),score(:,3)*std(STN_vx(:,3)),...
    'VariableNames',{'X','Y','Z','PC1','PC2','PC3'});

figure('renderer','painters','position',[300 300 1600 500])
tiledlayout(1,4)
nexttile
patch('faces',STN{2}.faces, 'vertices',STN{2}.vertices,"FaceColor",[0.9290 0.6940 0.1250],"edgecolor","none","FaceAlpha",0.4)
patch('faces',STN{3}.faces, 'vertices',STN{3}.vertices,"FaceColor",[0 0.4470 0.7410],"edgecolor","none","FaceAlpha",0.4)
patch('faces',STN{4}.faces, 'vertices',STN{4}.vertices,"FaceColor",[1 1 0],"edgecolor","none","FaceAlpha",0.4)
hold on
for pc_i = 1:3
    plot3(mean(STN_vx(:,1)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(1,pc_i),1,2).*[-1 1], mean(STN_vx(:,2)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(2,pc_i),1,2).*[-1 1], mean(STN_vx(:,3)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(3,pc_i),1,2).*[-1 1],'linewidth',4,'color','r');
    %plot3(mean(ClusterE_MNI_valid(:,1)) + coeff(1,pc_i) + param_mapping, mean(ClusterE_MNI_valid(:,2))+ coeff(2,pc_i) + param_mapping, mean(ClusterE_MNI_valid(:,3)) + coeff(3,pc_i) + param_mapping,'linewidth',2,'color',color_pc(pc_i,:))
end
set(gca,"view",1.0e+02 *[ -1.394536101759164   0.089355260649410])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
box off

for pc_i = 1:3
    nexttile
    s = scatter3(STN_PCtbl.X,STN_PCtbl.Y,STN_PCtbl.Z,15,STN_PCtbl.(['PC',num2str(pc_i)]),'filled');
    %s.SizeData = PC_tbl.nPoints;
    colormap(linspecer)
    view(-60,45)
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    cb = colorbar;
    ylabel(cb,sprintf('PC%d score [a.u.]',pc_i))
    box off
    set(gca,"view",1.0e+02 *[ -1.394536101759164   0.089355260649410])
    grid off
end


ClusterLoc_PC = table();
ClusterLoc_PC.X_MNI = ClusterS_MNI(:,1);
ClusterLoc_PC.Y_MNI = ClusterS_MNI(:,2);
ClusterLoc_PC.Z_MNI = ClusterS_MNI(:,3);

tempPC = nan(height(ClusterLoc_PC),3);
for ii = 1 : height(ClusterLoc_PC)
    tt = ClusterLoc_PC{ii,{'X_MNI','Y_MNI','Z_MNI'}};
    tempPC(ii,:) = (tt - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);
end
ClusterLoc_PC.PC1 = tempPC(:,1);
ClusterLoc_PC.PC2 = tempPC(:,2);
ClusterLoc_PC.PC3 = tempPC(:,3);


% ##########################3
%   this is wrong!!!!!!!
% idx_PC2clust= dsearchn(STN_vx,delaunayn(STN_vx),ClusterLoc_PC{:,{'X_MNI','Y_MNI','Z_MNI'}});
% 
% ClusterLoc_PC.PC1 = STN_PCtbl.PC1(idx_PC2clust);
% ClusterLoc_PC.PC2 = STN_PCtbl.PC2(idx_PC2clust);
% ClusterLoc_PC.PC3 = STN_PCtbl.PC3(idx_PC2clust);
% #############################


%% PC coordinates plot

CLUSTER_PC = [ClusterLoc_PC.PC1 ClusterLoc_PC.PC2 ClusterLoc_PC.PC3];
CoordBands = cell(1,6);
for fi = 2:6
    CoordBands{fi} = CLUSTER_PC(ClusterCentroidBand == fi,:);
end

pairs_test = nchoosek(2:6,2);
p_test = cell(1,3);
diff_test = cell(1,3);
for xi = 1 : 3
    p_test{xi} = nan(size(pairs_test,1),1);
    diff_test{xi} =  nan(size(pairs_test,1),1);
    for pi = 1 :  size(pairs_test,1)
        [p_test{xi}(pi), diff_test{xi}(pi)]  = permutation_diffTest_indsamples(CoordBands{pairs_test(pi,1)}(:,xi), CoordBands{pairs_test(pi,2)}(:,xi), 500);
    end
end


figure('position',[200 200 900 600])
tiledlayout(5,3)
for fi = 2: 6
    nexttile
    raincloud_plot(CoordBands{fi}(:,1), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,1)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,1))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

        nexttile
    raincloud_plot(CoordBands{fi}(:,2), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,2)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

        nexttile
    raincloud_plot(CoordBands{fi}(:,3), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,3)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,3))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

end

%%




PairLoc_PC = table();
PairLoc_PC.X_MNI = DB.Pairs.Location{:,{'S_MNI_X'}};
PairLoc_PC.Y_MNI = DB.Pairs.Location{:,{'S_MNI_Y'}};
PairLoc_PC.Z_MNI = DB.Pairs.Location{:,{'S_MNI_Z'}};


tempPC = nan(height(PairLoc_PC),3);
for ii = 1 : height(PairLoc_PC)
    tt = PairLoc_PC{ii,{'X_MNI','Y_MNI','Z_MNI'}};
    tempPC(ii,:) = (tt - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);
end
PairLoc_PC.PC1 = tempPC(:,1);
PairLoc_PC.PC2 = tempPC(:,2);
PairLoc_PC.PC3 = tempPC(:,3);



% #####################################3
% this is wrong!!!!!!!!
% idx_PC2pair = dsearchn(STN_vx,delaunayn(STN_vx),DB.Pairs.Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}});
% PairLoc_PC.PC1 = STN_PCtbl.PC1(idx_PC2pair);
% PairLoc_PC.PC2 = STN_PCtbl.PC2(idx_PC2pair);
% PairLoc_PC.PC3 = STN_PCtbl.PC3(idx_PC2pair);
% #####################################3



clust_prob_SPCs_bands = cell(1,3); %3 is number PC, 6 is number freq bands

clust_prob_SPCs_bands_pcaxis = cell(1,3);
for pc_i = 1:3
    clust_prob_SPCs_bands_pcaxis{pc_i} = unique(ClusterLoc_PC.(['PC',num2str(pc_i)]));
end

smooth_f = 0.5%.25;
for pc_i = 1:3
    clust_prob_SPCs_bands{pc_i} = nan(6,numel(clust_prob_SPCs_bands_pcaxis{pc_i}),numel(T));
    for fi = 2:6
        fprintf('Analyzing pc-%d for freq %d \n', pc_i, fi)
        for loc_i = 1: numel(clust_prob_SPCs_bands_pcaxis{pc_i})
            for ti = 1 : numel(T)
                try
                    clust_prob_SPCs_bands{pc_i}(fi,loc_i,ti) = sum(cluster_prob(ClusterCentroidBand == fi & (ClusterLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & ClusterLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)),ti))/sum((PairLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & PairLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)));
                catch
                    clust_prob_SPCs_bands{pc_i}(fi,loc_i,ti) = nan;
                end
            end
        end
    end
end


figure('renderer','painters','position',[300 10 1400 1600])
tiledlayout(5,6)
for pc_i = 1:3
    intpoints = linspace(min(clust_prob_SPCs_bands_pcaxis{pc_i}), max(clust_prob_SPCs_bands_pcaxis{pc_i}), numel(clust_prob_SPCs_bands_pcaxis{pc_i}));
    for fi = 2:6
        
        int = interp1(clust_prob_SPCs_bands_pcaxis{pc_i}',100*squeeze(clust_prob_SPCs_bands{pc_i}(fi,:,:)), intpoints);
        % perform plot
        nexttile(1+2*(pc_i-1)+6*(fi-2),[1 2])
        imagesc(T,intpoints,int)
        hold on

        colormap(linspecer)
     
        set(gca,'TickLabelInterpreter','none')
        cb = colorbar;
        ylabel(cb,'Overlap Cluster [%]')
        box off
        for evt = [2 3 4 5]'
            plot(EvtTimes(evt)' *[1,1], [min(clust_prob_SPCs_bands_pcaxis{pc_i}) max(clust_prob_SPCs_bands_pcaxis{pc_i})], '--', 'Color', 'w', 'LineWidth', 1)
        end
        set(gca, 'XTick', EvtTimes([2 3 4 5]))
        xLabeling = EvtTypes(2:end-1);
        set(gca, 'XTickLabel', xLabeling)
        xtickangle(30)
        set(gca,'ColorScale','log')
        caxis([0 1.05])
        %title(bands.symbol(fi))
        %
        set(gca,'YDir','normal')
        switch pc_i
            case 1
                ylabel('PC1 [mm] motor. <-> assoc.')
            case 2
                ylabel('PC2 [mm] ventral <-> dorsal')
            case 3
                ylabel('PC3 [mm] lateral <-> medial')
        end
    end
end


figure('renderer','painters','position',[300 10 1400 1600])
tiledlayout(5,6)
for pc_i = 1:3
    intpoints = linspace(min(clust_prob_SPCs_bands_pcaxis{pc_i}), max(clust_prob_SPCs_bands_pcaxis{pc_i}), numel(clust_prob_SPCs_bands_pcaxis{pc_i}));
    for fi = 2:6
        
        int = interp1(clust_prob_SPCs_bands_pcaxis{pc_i}',100*squeeze(clust_prob_SPCs_bands{pc_i}(fi,:,:)), intpoints);
        % perform plot
        nexttile(1+2*(pc_i-1)+6*(fi-2),[1 2])
        imagesc(T,intpoints,int)
        hold on

        colormap(linspecer)
     
        set(gca,'TickLabelInterpreter','none')
        cb = colorbar;
        ylabel(cb,'Overlap Cluster [%]')
        box off
        for evt = [2 3 4 5]'
            plot(EvtTimes(evt)' *[1,1], [min(clust_prob_SPCs_bands_pcaxis{pc_i}) max(clust_prob_SPCs_bands_pcaxis{pc_i})], '--', 'Color', 'w', 'LineWidth', 1)
        end
        set(gca, 'XTick', EvtTimes([2 3 4 5]))
        xLabeling = EvtTypes(2:end-1);
        set(gca, 'XTickLabel', xLabeling)
        xtickangle(30)
        caxis([0 5])
        %title(bands.symbol(fi))
        %
        set(gca,'YDir','normal')
        switch pc_i
            case 1
                ylabel('PC1 [mm] motor. <-> assoc.')
            case 2
                ylabel('PC2 [mm] ventral <-> dorsal')
            case 3
                ylabel('PC3 [mm] lateral <-> medial')
        end
    end
end

%% Repeat the analysis taking itno account ROIs

smooth_f = 0.2;%.25;
for pc_i = 1:3
    clust_prob_SPCs_bands{pc_i} = nan(6,numel(ROI_AtlasLabels), numel(clust_prob_SPCs_bands_pcaxis{pc_i}),numel(T));
    for fi = 2:4
        fprintf('Analyzing pc-%d for freq %d \n', pc_i, fi)
        for loc_i = 1: numel(clust_prob_SPCs_bands_pcaxis{pc_i})
                 if mod(ti,200) == 0
                    fprintf('Analyzing %d-%d (%1.2f) \n', ti, numel(T),ti/numel(T));
                end
            for ti = 1 : numel(T)
           
                for roi_i = [1 3 4 5]
                    try
                        clusts = sum(cluster_prob(contains(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i})' & ClusterCentroidBand == fi & (ClusterLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & ClusterLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)),ti));
                        pairs = sum(contains(DB.Pairs.Location.E_atlas_label_Destrieux,ROI_AtlasLabels{roi_i})' & (PairLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & PairLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)));
                        clust_prob_SPCs_bands{pc_i}(fi,roi_i,loc_i,ti) = clusts/pairs;
                    catch
                        clust_prob_SPCs_bands{pc_i}(fi,roi_i,loc_i,ti) = nan;
                    end
                end
            end
        end
    end
end





%% plot marginal distribution

figure('position',[100 100 900 900])
tiledlayout(5,6)
for fi = 2:6
    nexttile
    plot(T,squeeze(mean(clust_prob_SPCs_bands{1}(fi,:,:))))
    xlabel('T')
    ylabel('PPC1')
    ylim([0 0.013])
    nexttile
    plot(T,squeeze(mean(clust_prob_SPCs_bands{2}(fi,:,:))))
    xlabel('T')
    ylabel('PPC2') 
        ylim([0 0.013])

    nexttile
    plot(T,squeeze(mean(clust_prob_SPCs_bands{3}(fi,:,:))))
    xlabel('T')
    ylabel('PPC3')  
        ylim([0 0.013])

    nexttile
    plot(clust_prob_SPCs_bands_pcaxis{1},squeeze(mean(clust_prob_SPCs_bands{1}(fi,:,:),3)))
    xlabel('PPC1')
    ylabel('SPC') 
        ylim([0 0.026])

    nexttile
    plot(clust_prob_SPCs_bands_pcaxis{2},squeeze(mean(clust_prob_SPCs_bands{2}(fi,:,:),3)))
    xlabel('PPC2')
    ylabel('SPC') 
            ylim([0 0.026])
    nexttile
    plot(clust_prob_SPCs_bands_pcaxis{3},squeeze(mean(clust_prob_SPCs_bands{3}(fi,:,:),3)))
    xlabel('PPC3')
    ylabel('SPC') 
            ylim([0 0.026])

%         nexttile
%     plot(clust_prob_SPCs_bands_pcaxis{3},squeeze(mean(clust_prob_SPCs_bands{3}(fi,:,:),3)))
%     xlabel('PPC1')
%     ylabel('SPC') 
end




%%

figure('renderer','painters','position',[300 300 600 400])
for fi = 2:6
    nexttile
    %donut(EROIs_theta_n',EROIs_theta_unique)
    act_PC1 = squeeze(clust_prob_SPCs_bands{1}(fi,:,:)*100);
    act_PC2 = squeeze(clust_prob_SPCs_bands{2}(fi,:,:)*100);
    act_PC3 = squeeze(clust_prob_SPCs_bands{3}(fi,:,:)*100);
    plot(T,mean(act_PC1,'omitnan'),'linewidth',1.8,'color','k','linestyle','-')
    hold on
    plot(T,mean(act_PC2,'omitnan'),'linewidth',1.8,'color','k','linestyle','--')
    plot(T,mean(act_PC3,'omitnan'),'linewidth',1.8,'color','k','linestyle','-.')

    for evt = [2 3 4 5]'
        plot(EvtTimes(evt)' *[1,1], [0 max([mean(act_PC1,'omitnan') mean(act_PC2,'omitnan')])], '--', 'Color', 'r', 'LineWidth', 1)
    end
    set(gca, 'XTick', EvtTimes([2 3 4 5]))
    xLabeling = EvtTypes(2:end-1);
    set(gca, 'XTickLabel', xLabeling)
    xtickangle(30)
    ylabel('Overlap Cluster [%]')
    box off
    grid off
    legend({'PC1','PC2','PC3'},'location','northeast')
  
    ylim([0 8])
end

% check focality
figure('position',[200 200 1600 400])
tiledlayout(1,5)
for fi = 2:6
    nexttile
    coords = ClusterS_MNI(ClusterCentroidBand == fi,:);
    true_dist_clust = mean(vecnorm(coords - mean(coords,'omitnan'),2,2),'omitnan'); % mm
    perm_dist_clust = nan(1,nperms);
    for perm_i = 1 : nperms
        perm_loc = ClusterS_MNI(randi(height(ClusterS_MNI),1,size(coords,1)),:);
        perm_dist_clust(perm_i) = mean(vecnorm(perm_loc - mean(perm_loc,'omitnan'),2,2),'omitnan'); % mm
    end
    histogram(perm_dist_clust,round(sqrt(nperms)),'FaceColor',[.8 ,.8 ,.8])
    hold on
    box off
    xlabel('Distance to centre [mm]')
    ylabel(' # observations ')
    xlim([0.8 3])
    if true_dist_clust > min(perm_dist_clust)
        p_val = (sum((perm_dist_clust - mean(perm_dist_clust)) <= (true_dist_clust - mean(perm_dist_clust)))+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
    else
        p_val = 1/nperms;
    end
    xline(true_dist_clust,'label',{sprintf('p = %1.3f',p_val)},'Color','r','linewidth',1.4,'LabelOrientation','horizontal','FontSize',15)
end



%% check contribution different ROI
% theta
% clust 1: PC1: [-4.5 -3.5] (-4.09 peak) s; during speech production
%          PC2: [ 0 1] (2.03 peak) s; 
%          PC3: [0.8 1.3];
% alpha
% clust 1: PC1: [-4.5 -3.5]
%          PC2: [0 1]
%          PC3: [0.8 1.3];
% clust 2: PC1: [-2 -1.5];
%          PC2: [-3.5 -2.5]
%          PC3: [-1.3 -0.8]
% beta
% clust 1: PC1: [

PC1_int = [-4.5 -3.5];
PC2_int = [0 1];
PC3_int = [0.8 1.3];

% theta_map_PC1 = 100*squeeze(clust_prob_SPCs_bands{1}(fi,:,:));
% theta_map_PC2 = 100*squeeze(clust_prob_SPCs_bands{2}(fi,:,:));

% clust 1: PC1: [-4.25 -3.75] (-4.09 peak) s; during speech production
%          PC2: [1.75 2.25] (2.03 peak) s; 
idx_PC1Cluster = ClusterLoc_PC.PC1 >= PC1_int(1) & ClusterLoc_PC.PC1 <= PC1_int(2);
idx_PC2Cluster = ClusterLoc_PC.PC2 >= PC2_int(1) & ClusterLoc_PC.PC2 <= PC2_int(2);
idx_PC3Cluster = ClusterLoc_PC.PC3 >= PC3_int(1) & ClusterLoc_PC.PC3 <= PC3_int(2);


idx_PC1Pair = PairLoc_PC.PC1 >= PC1_int(1) & PairLoc_PC.PC1 <= PC1_int(2);
idx_PC2Pair = PairLoc_PC.PC2 >= PC2_int(1) & PairLoc_PC.PC2 <= PC2_int(2);
idx_PC3Pair = PairLoc_PC.PC3 >= PC3_int(1) & PairLoc_PC.PC3 <= PC3_int(2);


% theta peak
EROI_theta_peak_cluster = ClusterEAtlasLabel(ClusterCentroidBand' == fi & idx_PC1Cluster & idx_PC2Cluster & idx_PC3Cluster,:);
EROI_theta_peak_pair = ECoG_atlas(idx_PC1Pair & idx_PC2Pair & idx_PC3Pair,:);

EROI_theta_peak_density = cellfun(@(x) sum(strcmpi(EROI_theta_peak_cluster,x))/sum(strcmpi(EROI_theta_peak_pair,x))*100,ROI_AtlasLabels);
EROI_theta_peak_density(isnan(EROI_theta_peak_density)) = -5;
[~, idx_sort] = sort(EROI_theta_peak_density,'descend');

figure
nexttile
bar(1: numel(ROI_AtlasLabels),EROI_theta_peak_density(idx_sort))
xticks(1: numel(ROI_AtlasLabels))
xticklabels(ROI_AtlasLabels(idx_sort))
ylabel('t-SPC density')
box off

% first three are PostCG, STG and SMG

idx_roi = [5 4 3];
PairLoc_PC = struct2table(PairLoc_PC);
ThetaROI_SPCDensity = cell(1,numel(idx_roi));
for roi_i=  1 : numel(idx_roi)
    tmp = PairLoc_PC{strcmpi(ECoG_atlas,ROI_AtlasLabels{idx_roi(roi_i)}),{'X_MNI','Y_MNI','Z_MNI'}};
    [tmp_pairs,~,ix_pairs] = unique(tmp,'rows');
    % E_AtlasLabel_pairs_val_unique = E_AtlasLabel_pairs_val(ia_pairs);
    tmp_pairs_n = accumarray(ix_pairs,1)';
    tmp_clust_n = nan(size(tmp_pairs_n));
    clust_tmp = ClusterS_MNI(ClusterCentroidBand' == fi & strcmpi(ClusterEAtlasLabel, ROI_AtlasLabels{idx_roi(roi_i)}),:);
    for pair_i = 1  :numel(tmp_pairs_n)
        if tmp_pairs_n(pair_i) > 0
            tmp_clust_n(pair_i) = sum(ismember(clust_tmp,tmp_pairs(pair_i,:),'rows'));
        end
    end
    ThetaROI_SPCDensity{roi_i} = struct();
    ThetaROI_SPCDensity{roi_i}.loc = tmp_pairs;
    ThetaROI_SPCDensity{roi_i}.density = tmp_clust_n./tmp_pairs_n;

    idx_ThetaROIpair = dsearchn(STN_vx,delaunayn(STN_vx),ThetaROI_SPCDensity{roi_i}.loc);
    ThetaROI_SPCDensity{roi_i}.pcloc = [STN_PCtbl.PC1(idx_ThetaROIpair) STN_PCtbl.PC2(idx_ThetaROIpair) STN_PCtbl.PC3(idx_ThetaROIpair)];
end
    
ppk = nan(3,3);
ppm = nan(3,3);
ppkpc = nan(3,3);
ppmpc = nan(3,3);
for pp = 1:3
   [~,iii] = max(ThetaROI_SPCDensity{pp}.density);
   ppk(pp,:) = ThetaROI_SPCDensity{roi_i}.loc(iii,:);
   ppm(pp,:) = sum(ThetaROI_SPCDensity{pp}.density'.*ThetaROI_SPCDensity{pp}.loc,1)/sum(ThetaROI_SPCDensity{pp}.density);
   ppkpc(pp,:) = ThetaROI_SPCDensity{roi_i}.pcloc(iii,:);
   ppmpc(pp,:) = sum(ThetaROI_SPCDensity{pp}.density'.*ThetaROI_SPCDensity{pp}.pcloc,1)/sum(ThetaROI_SPCDensity{pp}.density);
end



%% run over all frequnecy and roi (if possible)
SPCROI_summary = struct();

% create nodz
% nodz_theta = cell2table(cell(0,5),'VariableNames',{'X','Y','Z','Color','Size'});
% nodz_alpha = cell2table(cell(0,5),'VariableNames',{'X','Y','Z','Color','Size'});
% nodz_beta = cell2table(cell(0,5),'VariableNames',{'X','Y','Z','Color','Size'});
nodz_theta = [];
nodz_alpha = [];
nodz_beta = [];


for roi_i = 1 : numel(ROI_AtlasLabels)
    tmp = PairLoc_PC{strcmpi(ECoG_atlas,ROI_AtlasLabels{roi_i}),{'X_MNI','Y_MNI','Z_MNI'}};
    [tmp_pairs,~,ix_pairs] = unique(tmp,'rows');
    % E_AtlasLabel_pairs_val_unique = E_AtlasLabel_pairs_val(ia_pairs);
    tmp_pairs_n = accumarray(ix_pairs,1)';
    SPCROI_summary(roi_i).roi = ROI_AtlasLabels{roi_i};
    SPCROI_summary(roi_i).X = tmp_pairs(:,1);
    SPCROI_summary(roi_i).Y = tmp_pairs(:,2);
    SPCROI_summary(roi_i).Z = tmp_pairs(:,3);
    
    % get pc coordinates
    PC = (((tmp_pairs - mean(STN_vx))./std(STN_vx))*coeff).*std(STN_vx);
    SPCROI_summary(roi_i).PC1 = PC(:,1);
    SPCROI_summary(roi_i).PC2 = PC(:,2);
    SPCROI_summary(roi_i).PC3 = PC(:,3);
    
    SPCROI_summary(roi_i).cog = nan(5,3);
    SPCROI_summary(roi_i).PCcog = nan(5,3);
    SPCROI_summary(roi_i).peak = nan(5,3);
    SPCROI_summary(roi_i).PCpeak = nan(5,3);

    for fi = 2:6
        tmp_clust_n = nan(size(tmp_pairs_n));
        clust_tmp = ClusterS_MNI(ClusterCentroidBand' == fi & strcmpi(ClusterEAtlasLabel, ROI_AtlasLabels{roi_i}),:);
        for pair_i = 1  :numel(tmp_pairs_n)
            if tmp_pairs_n(pair_i) > 0
                tmp_clust_n(pair_i) = sum(ismember(clust_tmp,tmp_pairs(pair_i,:),'rows'));
            end
        end
        SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ','')) = tmp_clust_n./tmp_pairs_n;
        % metrics
        SPCROI_summary(roi_i).cog(fi-1,:) = sum(SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ',''))'.*tmp_pairs,1,'omitnan')/sum(SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ','')));
        SPCROI_summary(roi_i).PCcog(fi-1,:) = sum(SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ',''))'.*PC,1,'omitnan')/sum(SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ','')));
        [~,iii] = max(SPCROI_summary(roi_i).(strrep(char(bands.name(fi)),' ','')));
        SPCROI_summary(roi_i).peak(fi-1,:)  = tmp_pairs(iii,:);
        SPCROI_summary(roi_i).PCpeak(fi-1,:) = PC(iii,:);
    end
    
    nodz_theta = [nodz_theta; [SPCROI_summary(roi_i).X SPCROI_summary(roi_i).Y SPCROI_summary(roi_i).Z roi_i*ones(numel(SPCROI_summary(roi_i).X),1) SPCROI_summary(roi_i).theta']];
    nodz_alpha = [nodz_alpha; [SPCROI_summary(roi_i).X SPCROI_summary(roi_i).Y SPCROI_summary(roi_i).Z roi_i*ones(numel(SPCROI_summary(roi_i).X),1) SPCROI_summary(roi_i).alpha']];
    nodz_beta= [nodz_beta; [SPCROI_summary(roi_i).X SPCROI_summary(roi_i).Y SPCROI_summary(roi_i).Z roi_i*ones(numel(SPCROI_summary(roi_i).X),1) SPCROI_summary(roi_i).beta']];

end

% nodz_theta(:,5) = (log10(nodz_theta(:,5) + eps) - min(log10(nodz_theta(:,5)+ eps)))/10;
% nodz_alpha(:,5) = (log10(nodz_alpha(:,5)+ eps) - min(log10(nodz_alpha(:,5)+ eps)))/10;
% nodz_beta(:,5) = (log10(nodz_beta(:,5)+ eps) - min(log10(nodz_beta(:,5)+ eps)))/10;

nodz_theta(nodz_theta(:,5) == 0 ,:) = [];
nodz_alpha(nodz_alpha(:,5) == 0 ,:) = [];
nodz_beta(nodz_beta(:,5) == 0 ,:) = [];

nodz_theta(:,5) = log10(nodz_theta(:,5)) - min(log10(nodz_theta(:,5)));
nodz_alpha(:,5) = log10(nodz_alpha(:,5)) - min(log10(nodz_alpha(:,5)));
nodz_beta(:,5) = log10(nodz_beta(:,5)) - min(log10(nodz_beta(:,5)));

writetable(array2table(nodz_theta,'VariableNames',{'X','Y','Z','Color','Size'}),'theta_cortical_projections.txt');
writetable(array2table(nodz_alpha,'VariableNames',{'X','Y','Z','Color','Size'}),'alpha_cortical_projections.txt');
writetable(array2table(nodz_beta,'VariableNames',{'X','Y','Z','Color','Size'}),'beta_cortical_projections.txt');

copyfile('theta_cortical_projections.txt','theta_cortical_projections.node')
copyfile('alpha_cortical_projections.txt','alpha_cortical_projections.node')
copyfile('beta_cortical_projections.txt','beta_cortical_projections.node')


% cretae nifti


Theta_cog = nan(numel(ROI_AtlasLabels),3);
Theta_peak = nan(numel(ROI_AtlasLabels),3);
Theta_PCcog = nan(numel(ROI_AtlasLabels),3);
Theta_PCpeak = nan(numel(ROI_AtlasLabels),3);
Alpha_cog = nan(numel(ROI_AtlasLabels),3);
Alpha_peak = nan(numel(ROI_AtlasLabels),3);
Alpha_PCcog = nan(numel(ROI_AtlasLabels),3);
Alpha_PCpeak = nan(numel(ROI_AtlasLabels),3);
Beta_cog = nan(numel(ROI_AtlasLabels),3);
Beta_peak = nan(numel(ROI_AtlasLabels),3);
Beta_PCcog = nan(numel(ROI_AtlasLabels),3);
Beta_PCpeak = nan(numel(ROI_AtlasLabels),3);

for roi_i = 1 : numel(ROI_AtlasLabels)
    Theta_cog(roi_i,:) = SPCROI_summary(roi_i).cog(1,:);
    Alpha_cog(roi_i,:) = SPCROI_summary(roi_i).cog(2,:);
    Beta_cog(roi_i,:) = SPCROI_summary(roi_i).cog(3,:);

    Theta_PCcog(roi_i,:) = SPCROI_summary(roi_i).PCcog(1,:);
    Alpha_PCcog(roi_i,:) = SPCROI_summary(roi_i).PCcog(2,:);
    Beta_PCcog(roi_i,:) = SPCROI_summary(roi_i).PCcog(3,:);

    Theta_peak(roi_i,:) = SPCROI_summary(roi_i).peak(1,:);
    Alpha_peak(roi_i,:) = SPCROI_summary(roi_i).peak(2,:);
    Beta_peak(roi_i,:) = SPCROI_summary(roi_i).peak(3,:);

    Theta_PCpeak(roi_i,:) = SPCROI_summary(roi_i).PCpeak(1,:);
    Alpha_PCpeak(roi_i,:) = SPCROI_summary(roi_i).PCpeak(2,:);
    Beta_PCpeak(roi_i,:) = SPCROI_summary(roi_i).PCpeak(3,:);
end

ROI_AtlasLabels_new = {'PreCG','MFG','SMG','STG','PostCG','SCG','IFG'};
theta_idx = [3 4 5];
alpha_idx = [2 3 4 5];
beta_idx = [1 3 5 6];

figure
tiledlayout(3,4)
nexttile
scatter3(Theta_cog(theta_idx,1),Theta_cog(theta_idx,2),Theta_cog(theta_idx,3),'filled')
hold on
text(Theta_cog(theta_idx,1), Theta_cog(theta_idx,2), Theta_cog(theta_idx,3), ROI_AtlasLabels_new(theta_idx))

nexttile
scatter3(Theta_PCcog(theta_idx,1),Theta_PCcog(theta_idx,2),Theta_PCcog(theta_idx,3),'filled')
hold on
text(Theta_PCcog(theta_idx,1), Theta_PCcog(theta_idx,2), Theta_PCcog(theta_idx,3), ROI_AtlasLabels_new(theta_idx))

nexttile
scatter3(Theta_peak(theta_idx,1),Theta_peak(theta_idx,2),Theta_peak(theta_idx,3),'filled')
hold on
text(Theta_peak(theta_idx,1), Theta_peak(theta_idx,2), Theta_peak(theta_idx,3), ROI_AtlasLabels_new(theta_idx))

nexttile
scatter3(Theta_PCpeak(theta_idx,1),Theta_PCpeak(theta_idx,2),Theta_PCpeak(theta_idx,3),'filled')
hold on
text(Theta_PCpeak(theta_idx,1), Theta_PCpeak(theta_idx,2), Theta_PCpeak(theta_idx,3), ROI_AtlasLabels_new(theta_idx))

nexttile
scatter3(Alpha_cog(alpha_idx,1),Alpha_cog(alpha_idx,2),Alpha_cog(alpha_idx,3),'filled')
hold on
text(Alpha_cog(alpha_idx,1), Alpha_cog(alpha_idx,2), Alpha_cog(alpha_idx,3), ROI_AtlasLabels_new(alpha_idx))

nexttile
scatter3(Alpha_PCcog(alpha_idx,1),Alpha_PCcog(alpha_idx,2),Alpha_PCcog(alpha_idx,3),'filled')
hold on
text(Alpha_PCcog(alpha_idx,1), Alpha_PCcog(alpha_idx,2), Alpha_PCcog(alpha_idx,3), ROI_AtlasLabels_new(alpha_idx))

nexttile
scatter3(Alpha_peak(alpha_idx,1),Alpha_peak(alpha_idx,2),Alpha_peak(alpha_idx,3),'filled')
hold on
text(Alpha_peak(alpha_idx,1), Alpha_peak(alpha_idx,2), Alpha_peak(alpha_idx,3), ROI_AtlasLabels_new(alpha_idx))

nexttile
scatter3(Alpha_PCpeak(alpha_idx,1),Alpha_PCpeak(alpha_idx,2),Alpha_PCpeak(alpha_idx,3),'filled')
hold on
text(Alpha_PCpeak(alpha_idx,1), Alpha_PCpeak(alpha_idx,2), Alpha_PCpeak(alpha_idx,3), ROI_AtlasLabels_new(alpha_idx))

nexttile
scatter3(Beta_cog(beta_idx,1),Beta_cog(beta_idx,2),Beta_cog(beta_idx,3),'filled')
hold on
text(Beta_cog(beta_idx,1), Beta_cog(beta_idx,2), Beta_cog(beta_idx,3), ROI_AtlasLabels_new(beta_idx))

nexttile
scatter3(Beta_PCcog(beta_idx,1),Beta_PCcog(beta_idx,2),Beta_PCcog(beta_idx,3),'filled')
hold on
text(Beta_PCcog(beta_idx,1), Beta_PCcog(beta_idx,2), Beta_PCcog(beta_idx,3), ROI_AtlasLabels_new(beta_idx))

nexttile
scatter3(Beta_peak(beta_idx,1),Beta_peak(beta_idx,2),Beta_peak(beta_idx,3),'filled')
hold on
text(Beta_peak(beta_idx,1), Beta_peak(beta_idx,2), Beta_peak(beta_idx,3), ROI_AtlasLabels_new(beta_idx))

nexttile
scatter3(Beta_PCpeak(beta_idx,1),Beta_PCpeak(beta_idx,2),Beta_PCpeak(beta_idx,3),'filled')
hold on
text(Beta_PCpeak(beta_idx,1), Beta_PCpeak(beta_idx,2), Beta_PCpeak(beta_idx,3), ROI_AtlasLabels_new(beta_idx))

%%
% figure
% for roi_i = 1 : 3
% nexttile
% scatter3(ThetaROI_SPCDensity{roi_i}.loc(:,1), ThetaROI_SPCDensity{roi_i}.loc(:,2), ThetaROI_SPCDensity{roi_i}.loc(:,3),ThetaROI_SPCDensity{roi_i}.density+.01,'filled')
% end
clust_prob_thetaROI = cell(numel(idx_roi),3); 

clust_prob_thetaROI_pcaxis = cell(1,3);
for pc_i = 1:3
    clust_prob_thetaROI_pcaxis{pc_i} = unique(ClusterLoc_PC.(['PC',num2str(pc_i)]));
end

smooth_f = 0.2%.25;
for roi_i = 1 : numel(idx_roi)
for pc_i = 1:3
    clust_prob_thetaROI{roi_i,pc_i} = nan(6,numel(clust_prob_thetaROI_pcaxis{pc_i}),numel(T));
  
        fprintf('Analyzing pc-%d for roi %d \n', pc_i, roi_i)
        for loc_i = 1: numel(clust_prob_thetaROI_pcaxis{pc_i})
            for ti = 1 : numel(T)
                try
                    clust_prob_thetaROI{roi_i,pc_i}(fi,loc_i,ti) = sum(cluster_prob(ClusterCentroidBand == fi & (ClusterLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & ClusterLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)),ti))/sum((PairLoc_PC.(['PC',num2str(pc_i)])' <= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) + smooth_f) & PairLoc_PC.(['PC',num2str(pc_i)])' >= (clust_prob_SPCs_bands_pcaxis{pc_i}(loc_i) - smooth_f)));
                catch
                    clust_prob_thetaROI{roi_i,pc_i}(fi,loc_i,ti) = nan;
                end
            end
        end
end
end

%% target relationship
figure('position',[100 100 800 800])
tiledlayout(5,5)
r_all = nan(5,5);
p_all = nan(5,5);
for fi = 2:6
    coords = SPCDensity{:,{'X_MNI' ,'Y_MNI','Z_MNI'}};
    value = SPCDensity{:,3+fi};
    nexttile
    dist = vecnorm(coords - targets.Alpha_Horn2017,2,2);
    [r_all(fi,1),p_all(fi,1)] = corr(dist,value);
    scatter(dist, value)
    nexttile
    dist = vecnorm(coords - targets.Beta_Horn2017,2,2);
    scatter(dist, value)
    [r_all(fi,2),p_all(fi,2)] = corr(dist,value);
    nexttile
    dist = vecnorm(coords - targets.BetaCoh_Wijk2022,2,2);
    scatter(dist, value)
    [r_all(fi,5),p_all(fi,5)] = corr(dist,value);

    nexttile
    dist = vecnorm(coords - targets.DBSOpt_Caire2013,2,2);
    scatter(dist, value)
    [r_all(fi,4),p_all(fi,4)] = corr(dist,value);

    nexttile
    dist = vecnorm(coords - targets.ThetaAlphaCoh_Wijk2022,2,2);
    scatter(dist, value)
    [r_all(fi,5),p_all(fi,5)] = corr(dist,value);

end




%%

% ClusterTheta = mean(clust_prob(ClusterLoc_PC.PC1 <=  - 3.75 & ClusterLoc_PC.PC1 >=  -4.25 & ClusterLoc_PC.PC2 <=  2.25 & ClusterLoc_PC.PC2 >=  1.75,:))/size(clust_prob,1);
% 
% figure
% plot(T,ClusterTheta)
% figname = fullfile(PATHS.saveFigures,'group-level','pca_scores_STN');
% saveas(gcf,figname,'fig')
% exportgraphics(gcf,[figname,'.png'])
% exportgraphics(gcf,[figname,'.eps'])
% fh{1} = figure;
% for fi = 1:5
%     idx_ = SPCDensity.Color == fi;
%     hold on
%     scatter(SPCDensity.Y(idx_), SPCDensity.Z(idx_),SPCDensity.Size(idx_)*500, ...
%         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerFaceAlpha',.3)
%     scatter(ClusterS_MNI_Band_avg{fi}(2),ClusterS_MNI_Band_avg{fi}(3), 400, ...
%         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerFaceAlpha',1,'Marker','^','linewidth',3)
% end