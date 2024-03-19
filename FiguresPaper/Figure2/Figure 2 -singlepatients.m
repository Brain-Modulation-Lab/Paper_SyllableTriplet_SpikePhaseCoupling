clc
clear all
close all

%% Step 0: settings of the script
Figure_2_settings;

%% Step 1. gather all SPC (pairs & clusters) information
DB_all = get_SPCdata(SUBJECTS,'all', PATHS);

%% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_all.Pairs);

%% Step 3. Get PPC maps single patients characteristics
PPCmap_mods = bml_annot_read_tsv(fullfile(PATH_PROJECT,'patients','Singlepatients_PPCmaps.tsv'));

%% Step 4: extract maps for each group
Groups = {'Lowfreq','Beta','No'};
nGroups = numel(Groups);
SPCmap = struct();

for groups_i = 1 : nGroups
    subj_group = PPCmap_mods.ID(PPCmap_mods.(Groups{groups_i}) == 1);
    idx_group = find(contains(SUBJECTS,subj_group))';
    DB = get_SPCdata(SUBJECTS,idx_group, PATHS);
    %DB = get_SPCdata(SUBJECTS,[1:5 7:25], PATHS);
    % 3.1a Show PPCz on Speech
    T = [-2.5 2];
    SPCmap.(Groups{groups_i}) = get_SPCmap(DB.Pairs,'PPCz', T,Tres,cfg,'Speech','all');
    SPCmap.(Groups{groups_i}).subjects = subj_group;
    clear DB

    cfg_fh = struct();
    cfg_fh.Position = [300 300 700 500];
    cfg_fh.EvtTimes = EvtTimes.Speech;
    cfg_fh.EvtTypes = cfg.plot.EvtTypes;
    %cfg_fh.CLim = [-.05 .15];
    cfg_fh.CMap = linspecer;

    fh = plotter_SPCmap(SPCmap.(Groups{groups_i}), cfg_fh);
    figname = fullfile('W:\Users\MV1019\PhaseLocking\FiguresPaper\Figure2\patients',['PPCzmap_pairs-all_evt-speech_subj-', Groups{groups_i}]);
    saveFigures(fh, figname)
    close gcf
end

%% Step 4.5: get ROIs
% Use Destrieux atlas
atlas = 'E_atlas_label_Destrieux';
min_subj = 7;
min_pairs = 100;
flag_sort = true;
[ROI_AtlasLabels,ROIAtlasLabel_pairs_n]  = getROIs(DB_all.Pairs.Location.(atlas), cellstr(DB_all.Pairs.Location.subj_id), min_subj,min_pairs,flag_sort);

%% Step 5. Get STN and ECoG coordinates (& Atlas labels) for the patients in the three group
res_sphere= 2; % mm
for groups_i = 1 : nGroups
    subj_group = PPCmap_mods.ID(PPCmap_mods.(Groups{groups_i}) == 1);
    idx_group = find(contains(SUBJECTS,subj_group))';
    % get data
    DB = get_SPCdata(SUBJECTS,idx_group, PATHS);
    % get coordinates
    Coords = get_SPCcoords(DB.Pairs.Location);

    % create nifti and nodes ECOG files 
    NIFTINAME = fullfile(PATH_OUTPUT,'figures',TYPE_DB, ...
        ['Coverage_area-ECoG_res-',num2str(res_sphere), ...
        'mm_map-Npairs_cond-',Groups{groups_i},'.nii']);
    create_Locmap(Coords.ECoG, res_sphere, CTX,NIFTINAME);
    NodeMap = LocMap2NodeMap(NIFTINAME);

    % create nifti and nodes STN files 
    NIFTINAME = fullfile(PATH_OUTPUT,'figures',TYPE_DB, ...
        ['Coverage_area-STN_res-',num2str(res_sphere), ...
        'mm_map-Npairs_cond-',Groups{groups_i},'.nii']);
    create_Locmap(Coords.STN, res_sphere, STN_nii,NIFTINAME);
    NodeMap = LocMap2NodeMap(NIFTINAME);
end

