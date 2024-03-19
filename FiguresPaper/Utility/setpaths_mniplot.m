function [atlases, average_mni, BS1] = setpaths_mniplot(PATH_DISTAL,PATH_AVERAGE_MNI, PATH_CORTEX_MNI)

load(fullfile(PATH_DISTAL,'atlas_index.mat'),'atlases'); % manually load definition of DISTAL atlas.
average_mni = load(PATH_AVERAGE_MNI);
load(PATH_CORTEX_MNI,'BS1');
