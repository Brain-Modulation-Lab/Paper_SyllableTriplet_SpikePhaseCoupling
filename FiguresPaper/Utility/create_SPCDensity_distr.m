function fh = create_SPCDensity_distr(cfg,DB)

bands = cfg.bands;
nperms = cfg.nperms;
atlas = cfg.atlas;
min_subj = cfg.min_subj;
min_pairs = cfg.min_pairs;
flag_sort = cfg.flag_sort;
% load clusters
SPCDensity = readtable(cfg.map);
% load cortex
load(cfg.ref);

% figure()
% tiledlayout(5,3)
% for fi = 1:5
%     idx_ = SPCDensity.Color == fi;
% nexttile
% scatter(SPCDensity.X(idx_),SPCDensity.Size(idx_),15, ...
%         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerFaceAlpha',.8)
% nexttile
% scatter(SPCDensity.Y(idx_),SPCDensity.Size(idx_),15, ...
%         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerFaceAlpha',.8)
% nexttile
% scatter(SPCDensity.Z(idx_),SPCDensity.Size(idx_),15, ...
%         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         'MarkerFaceAlpha',.8)
% 
% end
    % find coordinates (Y for CS and Z for FS)

cfg.VarFields = {'ClusterPhase','ClusterCentroid','ClusterCentroidBand','ClusterE_MNI','ClusterEAtlasLabel','ClusterDur'};
[ClusterPhase,ClusterCentroid,ClusterCentroidBand, ClusterE_MNI,ClusterEAtlasLabel,ClusterDur] = get_SPCClustersProperties(cfg,DB.Clusters);


%%
fh = figure('Position',[200 200 700 900]);
tiledlayout(3,1)

E_XMNIPairs = DB.Pairs.Location.E_MNI_X;
E_XMNI_vec = min(E_XMNIPairs) : 0.5 : max(E_XMNIPairs);

E_XMNIDensity = nan(6,numel(E_XMNI_vec));
E_XMNIClusterDensity = nan(6,numel(E_XMNI_vec));
for fi = 2:6
    E_XMNIClusters = ClusterE_MNI(ClusterCentroidBand == fi,1);
    IdxXMNI_Pairs = rangesearch(E_XMNIPairs,E_XMNI_vec',2);
    IdxXMNI_Clusters = rangesearch(E_XMNIClusters,E_XMNI_vec',2);
    E_XMNIDensity(fi,:) = cellfun(@(x,y) numel(x)/numel(y)*100,IdxXMNI_Clusters,IdxXMNI_Pairs);
    E_XMNIDensity(fi,isinf(E_XMNIDensity(fi,:))) = nan;
    E_XMNIClusterDensity(fi,:) = cellfun(@(x,y) numel(x),IdxXMNI_Clusters);
    E_XMNIClusterDensity(fi,:) = E_XMNIClusterDensity(fi,:)/max(E_XMNIClusterDensity(fi,:),[],'omitnan');
end


nexttile
hold on
for fi = 2:6
    plot(E_XMNI_vec,E_XMNIDensity(fi,:)/max(E_XMNIDensity(fi,:),[],'omitnan'),'color',hex2rgb(bands.color(fi)),'linewidth',2)
end
%xlim([-55,45])
xlabel('X coords [mm]')


E_YMNIPairs = DB.Pairs.Location.E_MNI_Y;
E_YMNI_vec = min(E_YMNIPairs) : 0.5 : max(E_YMNIPairs);

E_YMNIDensity = nan(6,numel(E_YMNI_vec));
E_YMNIClusterDensity = nan(6,numel(E_YMNI_vec));
for fi = 2:6
    E_YMNIClusters = ClusterE_MNI(ClusterCentroidBand == fi,2);
    IdxYMNI_Pairs = rangesearch(E_YMNIPairs,E_YMNI_vec',2);
    IdxYMNI_Clusters = rangesearch(E_YMNIClusters,E_YMNI_vec',2);
    E_YMNIDensity(fi,:) = cellfun(@(x,y) numel(x)/numel(y)*100,IdxYMNI_Clusters,IdxYMNI_Pairs);
    E_YMNIDensity(fi,isinf(E_YMNIDensity(fi,:))) = nan;
    E_YMNIClusterDensity(fi,:) = cellfun(@(x,y) numel(x),IdxYMNI_Clusters);
    E_YMNIClusterDensity(fi,:) = E_YMNIClusterDensity(fi,:)/max(E_YMNIClusterDensity(fi,:),[],'omitnan');
end


nexttile
hold on
for fi = 2:6
    plot(E_YMNI_vec,E_YMNIDensity(fi,:)/max(E_YMNIDensity(fi,:),[],'omitnan'),'color',hex2rgb(bands.color(fi)),'linewidth',2)
end
xlim([-55,45])
xlabel('Y coords [mm]')


E_ZMNIPairs = DB.Pairs.Location.E_MNI_Z;
E_ZMNI_vec = min(E_ZMNIPairs) : 0.5 : max(E_ZMNIPairs);

E_ZMNIDensity = nan(6,numel(E_ZMNI_vec));
E_ZMNIClusterDensity = nan(6,numel(E_ZMNI_vec));
for fi = 2:6
    E_ZMNIClusters = ClusterE_MNI(ClusterCentroidBand == fi,3);
    IdxZMNI_Pairs = rangesearch(E_ZMNIPairs,E_ZMNI_vec',2);
    IdxZMNI_Clusters = rangesearch(E_ZMNIClusters,E_ZMNI_vec',2);
    E_ZMNIDensity(fi,:) = cellfun(@(x,y) numel(x)/numel(y)*100,IdxZMNI_Clusters,IdxZMNI_Pairs);
    E_ZMNIDensity(fi,isinf(E_ZMNIDensity(fi,:))) = nan;
    E_ZMNIClusterDensity(fi,:) = cellfun(@(x,y) numel(x),IdxZMNI_Clusters);
    E_ZMNIClusterDensity(fi,:) = E_ZMNIClusterDensity(fi,:)/max(E_ZMNIClusterDensity(fi,:),[],'omitnan');
end

nexttile
hold on
for fi = 2:6
    plot(E_ZMNI_vec,E_ZMNIDensity(fi,:)/max(E_ZMNIDensity(fi,:),[],'omitnan'),'color',hex2rgb(bands.color(fi)),'linewidth',2)
end
%xlim([-55,45])
xlabel('Z coords [mm]')

% for fi = 2:6
%     plot(E_YMNI_vec,E_YMNIClusterDensity(fi,:),'color',hex2rgb(bands.color(fi)))
% end
% 
% fh = figure('Position',[200 200 400 300]);
% hold on
% for fi = 2:6
%     scatter(SPCDensity.Y(SPCDensity.Color == (fi-1)),SPCDensity.Size(SPCDensity.Color == (fi-1)),'filled')
% end





