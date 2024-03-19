
% this script calculates spike sorting quality metrics
% 1. ISI violation (% < 3 ms)
% 2. mean firing rate
% 3. peak SNR
% 4. isolation distance
% 5. burst index

% need to define the neurons used in the analysus

PATH_SERVER = 'Z:\DBS';
if ispc
    PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking';
else
    PATH_PROJECT = '/Volumes/Nexus4/Users/MV1019/PhaseLocking';
end
PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
PATH_OUTPUT = fullfile(PATH_PROJECT,'Supplementary_Analysis');
PATH_DISTAL = 'Z:\Resources\DISTAL-Atlas';
PATH_NIFTI = fullfile(PATH_OUTPUT,'Nifti');

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
cfg = set_configs('default');



fmin  = 2;
fmax = 140;
nfreq = 60;

fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

TYPE_DB = 'default';
PATHS = struct();
PATHS.DataDir      = fullfile('Y:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_PROJECT,'groupanalyses','Output','Results_new',TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures');


bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;


% build information abotu unti types used in the analysis
% PairsLocation_MNI_all = [];
% toDo = [1:4 6:n_SUBJECTS];
% for subj_i = toDo
%
%     % TF = [];
%     SUBJECT = SUBJECTS{subj_i};
%       fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
%     PATH_clusters = fullfile(PATHS.saveMatFiles ,SUBJECT,'ClustersPLV.mat');
%
%
%     if isfile(PATH_clusters)
%         fprintf(' Stacking clusters results in %s \n', PATH_clusters)
%         load(PATH_clusters,'PairsLocation_MNI');
%         % put info about subjects id
%         PairsLocation_MNI.subj_id = repmat(SUBJECT, height(PairsLocation_MNI),1);
%         % stack location
%         PairsLocation_MNI_all = [PairsLocation_MNI_all; PairsLocation_MNI];
%     end
% end

%%

toDo = 1 : n_SUBJECTS;% 1 : n_SUBJECTS;
[nSign_pairs_all, nPairs_all]  = deal(nan(1,n_SUBJECTS)) ;
Clusters_all = [];
PPC_mat_all = [];
PPCz_mat_all = [];
PPCmedianperm_mat_all = [];
ES_mat_all = [];
Phase_mat_all = [];
PairsLocation_MNI_all = [];
PLVTimeEvts_all = [];
for subj_i = toDo %
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)

    PATH_clusters = fullfile(PATHS.saveMatFiles,SUBJECT,'ClustersPLV.mat');

    if isfile(PATH_clusters)
        fprintf(' Stacking clusters results in %s \n', PATH_clusters)
        load(PATH_clusters,'Clusters', 'nSign_pairs','n_pairs','PairsLocation_MNI', 'PPC_mat','PPCmedianperm_mat','PPCz_mat','ES_mat','PLVTimeEvts','Phase_mat');
        % put info about subjects id
        [Clusters.subj_id] = deal(SUBJECT);
        PairsLocation_MNI.subj_id = repmat(SUBJECT, height(PairsLocation_MNI),1);
        % stack Clusters
        Clusters_all = [Clusters_all Clusters(~isnan([Clusters.S_typeFRmod]))]; % eliminate nan firemod
        PairsLocation_MNI_all = [PairsLocation_MNI_all; PairsLocation_MNI];
        PPC_mat_all = [PPC_mat_all  PPC_mat];
        PPCz_mat_all = [PPCz_mat_all  PPCz_mat];
        PPCmedianperm_mat_all = [PPCmedianperm_mat_all  PPCmedianperm_mat];
        ES_mat_all = [ES_mat_all  ES_mat];
        Phase_mat_all = [Phase_mat_all Phase_mat];
        PLVTimeEvts_all = [PLVTimeEvts_all PLVTimeEvts];
        % grab information about significant pairs and n_pairs
        nSign_pairs_all(subj_i) = nSign_pairs;
        nPairs_all(subj_i) = n_pairs;

        fprintf(' Completed Stacking clusters results in %s \n', PATH_clusters)
        % clear variables for consistency
        clear Clusters nSign_pairs n_pairs PairsLocation_MNI PPC_mat PPCz_mat ES_mat PPCmedianperm_mat PLVTimeEvts Phase_mat
    else
        warning('Analysis still running: Clusters  is not available yet for Subject %s %d - %d \n', SUBJECT, subj_i, n_SUBJECTS)
    end
end

nClusters = numel(Clusters_all);
nonan = ~isnan(nSign_pairs_all); % to manage the lack of data for patient 5



unit_types = [-1 0 1 2];
unit_types_labels= {'Neg','None','Excit','Mixed'};
n_unit_types_labels = [40, 69 85 25];


%% extract unique units (211 neurons)
PairsLocation_MNI_all.subj_id = string(PairsLocation_MNI_all.subj_id);
uniqueUnits = unique([ PairsLocation_MNI_all.subj_id PairsLocation_MNI_all.S_channel],'rows');
Units = cellfun(@(x) uniqueUnits(contains(uniqueUnits(:,1),x),2), SUBJECTS,'uni',false);
% get sortnotes for each patients to understand where these neurons come
% from

for subj_i = 1 : n_SUBJECTS
    SUBJECT = SUBJECTS{subj_i};
    PATH_ANNOT = fullfile('Z:\DBS',SUBJECT,'Preprocessed Data\Sync\annot',[SUBJECT])

end