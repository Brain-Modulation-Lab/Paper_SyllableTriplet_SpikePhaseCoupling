function  Unit_Info = get_EncodingTypes(Loc,encoding_type,vargargin)

EncodingModel = load('model_fit_summary.mat');
EncodingTable = EncodingModel.rsq_table_select;


if istable(Loc)
    LocInfo = [Loc.subj_id Loc.S_channel];
elseif isstruct(Loc)
    if nargin  == 3
       Cluster_Subject = {Loc.subj_id}';
       Cluster_Schannel = {Loc.S_channel}';
    LocInfo = string([Cluster_Subject, Cluster_Schannel]);

    else
       cfg = [];
        cfg.VarFields = {'ClusterSubjects','ClusterSChannel'};
        [Cluster_Subject, Cluster_Schannel] = get_SPCClustersProperties(cfg,Loc);
       LocInfo = [Cluster_Subject, Cluster_Schannel];
 
    end
end



[~, crit_p, ~, ~]=fdr_bh(EncodingModel.rsq_data.p(:));
ttt = squeeze(EncodingModel.rsq_data.p <= 0.05)';
%ttt = EncodingTable(:,23:29);

switch encoding_type
    case 'consonant'
        idx_ = [1 2 3];
    case 'vowel'
        idx_ = [4 5];
    case 'position'
        idx_ = [6 7];
end

EncodingUnits_idx = find(any(ttt(:,idx_),2));
EncodingUnits_subj_id = EncodingTable.subject(EncodingUnits_idx);
S_channel = EncodingTable.channel(EncodingUnits_idx);
EncodingUnits_S_channel = cell(numel(EncodingUnits_idx),1);

Unit_Info = nan(height(LocInfo),1);

for ui = 1 : numel(EncodingUnits_idx)
    tmp = strsplit(S_channel{ui},'_');
    EncodingUnits_S_channel{ui} = [tmp{2},tmp{3}];
    Unit_Info(ismember(LocInfo,string([EncodingUnits_subj_id(ui) EncodingUnits_S_channel(ui)]),'rows')) = 1;
end







