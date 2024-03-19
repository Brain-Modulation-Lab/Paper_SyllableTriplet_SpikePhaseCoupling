if ismac
    PATH_SERVER = '/Volumes/Nexus/DBS';
    PATH_PROJECT = '/Volumes/Nexus4/Users/MV1019/PhaseLocking/FiguresPaper/Figures_4_5_spatial';
    PATH_RESOURCES = '/Volumes/Nexus1/Resources/MNI_Cortex_plotting';
    PATH_DISTAL = '/Volumes/Nexus/Resources/DISTAL-Atlas';
    PATH_UTILITIES =  '/Volumes/Nexus4/Users/MV1019/PhaseLocking/FiguresPaper/Utility';
elseif ispc
    PATH_SERVER = 'Z:\DBS';
    PATH_PROJECT = 'W:\Users\MV1019\PhaseLocking\FiguresPaper\Figures_4_5_spatial';
    PATH_RESOURCES = 'Z:\Resources\MNI_Cortex_plotting';
    PATH_DISTAL = 'Z:\Resources\DISTAL-Atlas';
    PATH_UTILITIES =  'W:\Users\MV1019\PhaseLocking\FiguresPaper\Utility';

end


PATH_OUTPUT = fullfile(PATH_PROJECT,'Output');

addpath(genpath(PATH_UTILITIES))

prettify_plots;
% ft_defaults;
% bml_defaults;

% laod subjects
if ispc
    SUBJECTS = readtable('W:\Users\MV1019\PhaseLocking\FiguresPaper\Subjects_list.txt');
elseif ismac
    SUBJECTS = readtable('/Volumes/Nexus4/Users/MV1019/PhaseLocking/FiguresPaper/Subjects_list.txt');
end

SUBJECTS = SUBJECTS.Subjects;
n_SUBJECTS = numel(SUBJECTS);

% Load time tol
if ispc
    TOL = bml_annot_read_tsv('W:\Users\MV1019\PhaseLocking\FiguresPaper\time-tolerance.tsv');
elseif ismac
    %TOL = bml_annot_read_tsv('/Volumes/Nexus4/Users/MV1019/PhaseLocking/FiguresPaper/time-tolerance.tsv');
end
% set paths

cfg = set_configs("default");
TYPE_DB = 'default';
fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;

PATHS = struct();
if ispc
    PATHS.DataDir      = fullfile('W:\Users\MV1019\PhaseLocking\','groupanalyses','Output','Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
elseif ismac
    PATHS.DataDir      = fullfile('/Volumes/Nexus4/Users/MV1019/PhaseLocking/','groupanalyses','Output','Results_new', TYPE_DB); % in this case they are the same because E and S are in the same folder of PLV
end
PATHS.saveMatFiles = fullfile(PATH_OUTPUT,'Results', TYPE_DB);
PATHS.saveFigures  = fullfile(PATH_OUTPUT,'Figures',TYPE_DB);


cfg.PATHS = PATHS;
% nifti settings
if exist('ea_space') == 2
    CTX = ea_load_nii(fullfile(ea_space,'labeling','AICHA (Joliot 2015).nii'));
    STN_nii = ea_load_nii(fullfile(ea_space,'atlases','DISTAL Minimal (Ewert 2017)','lh','STN.nii.gz'));

else
    if ispc
        load(fullfile('W:\Users\MV1019\PhaseLocking\groupanalyses\Output','Results_new','default','group-level','CTX.mat'))
    elseif ismac
        load(fullfile('/Volumes/Nexus4/Users/MV1019/PhaseLocking/groupanalyses/Output','Results_new','default','group-level','CTX.mat'))
    end
% nii = spm_vol('STN.nii');
% img = spm_read_vols(nii);
% volnum = numel(nii);
% STN_nii = nii(1); % only keep the header of the first vol for multi-vol image
% STN_nii.img = img; % set image data
% %STN_nii.voxsize = ea_detvoxsize(nii.mat); % set voxsize
% STN_nii.volnum = volnum; % set number of volumes

end



% parameters

bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;
Tres = 5E-3 * 10;
cfg.bands = bands;