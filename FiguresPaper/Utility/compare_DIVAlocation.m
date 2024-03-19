function fh = compare_DIVAlocation(cfg,DB)

bands = cfg.bands;
ref_coords = readtable(cfg.ref_coords);
SPCDensity = readtable(cfg.map);

ref_coords = ref_coords(~contains(ref_coords.Map,'respiratory') & ~contains(ref_coords.Map,'Initiation') & ~contains(ref_coords.Map,'Motivation'),:);
% Maps = unique(ref_coords.Map);
% Maps_coords = cell2mat(cellfun(@(x) mean(ref_coords{strcmpi(ref_coords.Map,x),{'X','Y','Z'}},1),Maps,'uni',false));

Maps = ref_coords.Map;
Maps_coords = ref_coords{:,{'X','Y','Z'}};
% get cluster coords
cfg.VarFields = {'ClusterCentroidBand','ClusterE_MNI'};
[ClusterCentroidBand,ClusterE_MNI] = get_SPCClustersProperties(cfg,DB.Clusters);

ClusterE_MNI_freq = arrayfun(@(x) ClusterE_MNI(ClusterCentroidBand == x,:),1:6,'uni',false);

Divadistances_mean = nan(6,numel(Maps));
Divadistances_sem = nan(6,numel(Maps));
Divadistances = nan(6,numel(Maps));

ClusterS_MNI_Band_peak = nan(5,3);

for fi = 1:5
    tmp = SPCDensity(SPCDensity.Color == fi,:);
    [~,idx_] = max(tmp.Size,[],'omitnan');
    ClusterS_MNI_Band_peak(fi,:) = tmp{idx_,{'X','Y','Z'}};
end

%Divadistances
for fi = 1 : 6
    if ~isempty(ClusterE_MNI_freq{fi})
        for map_i = 1 : numel(Maps)
            Divadistances_mean(fi,map_i) = mean(vecnorm(ClusterE_MNI_freq{fi} - Maps_coords(map_i,:),2,2),1,'omitnan');
           Divadistances_sem(fi,map_i) = std(vecnorm(ClusterE_MNI_freq{fi} - Maps_coords(map_i,:),2,2),[],1,'omitnan'); %/sqrt(numel(ClusterE_MNI_freq{fi}));
        end
    end
end
for fi = 2 : 6
    for map_i = 1 : numel(Maps)

        Divadistances(fi,map_i) = vecnorm(Maps_coords(map_i,:) - ClusterS_MNI_Band_peak(fi-1,:));
    end
end

fh = figure('position',[200 200 1600 600]);
imagesc(Divadistances(2:6,:))
colormap(linspecer)
colorbar
xticks(1:numel(Maps))
xticklabels(Maps)
yticks(1:5)
yticklabels(bands.symbol(2:6))
box off

% figure
% x = 1 : numel(Maps);
% for fi = 2:6
%     hold on
%   plot(x, Divadistances(fi,:) ...
%         ,'color',hex2rgb(bands.color(fi)),'linewidth',1.3);
% end
% xticks(x)
% xticklabels(Maps)
% legend(bands.symbol(2:6))
% box off


% figure
% x = 1 : numel(Maps);
% ph = [];
% for fi = 2:6
%     xShaded = [x fliplr(x)];
%     yShaded = [Divadistances_mean(fi,:) - Divadistances_sem(fi,:)   fliplr(Divadistances_mean(fi,:) + Divadistances_sem(fi,:))];
%     fill(xShaded, yShaded,hex2rgb(bands.color(fi)),'FaceAlpha',.3);
%     hold on
%     ph = [ ph plot(x, Divadistances_mean(fi,:) ...
%         ,'color',hex2rgb(bands.color(fi)),'linewidth',1.3)];
% end
% xticks(x)
% xticklabels(Maps)
% legend(ph,bands.symbol(2:6))
% box off

% % check locality dist(coords - peak) ~ SPC_density(coords)
% Rmaps =  nan(6,numel(Maps));
% Pmaps =  nan(6,numel(Maps));
% 
% for fi = 1:5
%     spc = SPCDensity{SPCDensity.Color == fi,{'Size'}};
%     coords = SPCDensity{SPCDensity.Color == fi,{'X','Y','Z'}};
%     for map_i = 1 : numel(Maps)
%         peak = Maps_coords(map_i,:);
% 
%         distance = sqrt((coords(:,1) - peak(1)).^2 + (coords(:,2) - peak(2)).^2 + (coords(:,3) - peak(3)).^2);
%         %     nexttile
%         %     scatter(-distance,spc,15, ...
%         %         'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         %         'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
%         %         'MarkerFaceAlpha',.3)
%         [Rmaps(fi+1,map_i),Pmaps(fi+1,map_i)] = corr(spc,-distance);
%         %     title(sprintf('R = %1.2f, p = %1.3f',r,p))
%         %     xlabel('Distance to Peak [mm]')
%         %     ylabel('SPC density')
% 
% %         if map_i == 1 && fi== 5
% %             figure
% %             scatter(distance,spc,15, ...
% %                 'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
% %                 'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
% %                 'MarkerFaceAlpha',.3)
% %             title(sprintf('R = %1.2f, p = %1.3f',Rmaps(fi+1,map_i),Pmaps(fi+1,map_i)))
% %             xlabel('Distance to Peak [mm]')
% %             ylabel('SPC density')
% %         end
% 
%         
%     end
% end
% figure('position',[200 200 1600 1300])
% tiledlayout(3,1)
% nexttile
% imagesc(Divadistances(2:6,:))
% colormap(linspecer)
% colorbar
% xticks(1:numel(Maps))
% xticklabels('')
% yticks(1:5)
% yticklabels(bands.symbol(2:6))
% box off
% nexttile
% imagesc(Rmaps(2:6,:))
% colormap(linspecer)
% colorbar
% xticks(1:numel(Maps))
% xticklabels('')
% yticks(1:5)
% yticklabels(bands.symbol(2:6))
% box off
% nexttile
% imagesc(Pmaps(2:6,:))
% colormap(linspecer)
% colorbar
% xticks(1:numel(Maps))
% xticklabels(Maps)
% yticks(1:5)
% yticklabels(bands.symbol(2:6))
% box off

% Divadistances_stat = Divadistances;
% Divadistances_stat(Pmaps >= 0.05 | Rmaps <= 0) = nan;
% figure('position',[200 200 1600 500])
% 
% 
% imagesc(Divadistances_stat(2:6,:))
% colormap(linspecer)
% colorbar
% xticks(1:numel(Maps))
% xticklabels(Maps)
% yticks(1:5)
% yticklabels(bands.symbol(2:6))
% box off