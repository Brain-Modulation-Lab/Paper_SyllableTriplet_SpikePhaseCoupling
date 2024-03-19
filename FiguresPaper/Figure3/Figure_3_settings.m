if ismac
    PATH_SERVER = '/Volumes/Nexus/DBS';
    PATH_PROJECT = '/Volumes/Nexus4/Users/MV1019/PhaseLocking';
    PATH_RESOURCES = '/Volumes/Nexus1/Resources/MNI_Cortex_plotting';
elseif ispc
    PATH_SERVER = 'Z:\DBS';
    PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking\FiguresPaper\Figure3';
    PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
end


PATH_OUTPUT = fullfile(PATH_PROJECT,'Output');

PATH_UTILITIES =  'W:\Users\MV1019\PhaseLocking\FiguresPaper\Utility';
addpath(genpath(PATH_UTILITIES))

prettify_plots;
ft_defaults;
bml_defaults;

% laod subjects
SUBJECTS = readtable('W:\Users\MV1019\PhaseLocking\FiguresPaper\Subjects_list.txt');
SUBJECTS = SUBJECTS.Subjects;
n_SUBJECTS = numel(SUBJECTS);

% Load time tol
TOL = bml_annot_read_tsv('W:\Users\MV1019\PhaseLocking\FiguresPaper\time-tolerance.tsv');

% set paths

cfg = set_configs("default");
TYPE_DB = 'default';
fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

PATHS = struct();
PATHS.DataDir      = fullfile('W:\Users\MV1019\PhaseLocking\','groupanalyses','Output','Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results', TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures',TYPE_DB);


cfg.PATHS = PATHS;
% nifti settings
if exist('ea_space') == 2
    CTX = ea_load_nii(fullfile(ea_space,'labeling','AICHA (Joliot 2015).nii'));
else
    load(fullfile('W:\Users\MV1019\PhaseLocking\groupanalyses\Output','Results_new','default','group-level','CTX.mat'))
end
STN_nii = ea_load_nii(fullfile(ea_space,'atlases','DISTAL Minimal (Ewert 2017)','lh','STN.nii.gz'));



% parameters

bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;
Tres = 5E-3 * 10;
cfg.bands = bands;