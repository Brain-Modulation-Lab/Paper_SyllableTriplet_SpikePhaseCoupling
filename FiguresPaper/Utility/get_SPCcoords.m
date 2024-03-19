function Coords = get_SPCcoords(Location,type)

if ~exist("type",'var')
    type = 'both';
end

switch type
    case 'both'
        Coords.STN = Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}};
        Coords.ECoG = Location{:,{'E_MNI_X','E_MNI_Y','E_MNI_Z'}};
        Coords.ECoG_atlas = Location{:,{;'E_atlas_label_Destrieux'}};   
    case 'stn'
        Coords.STN = Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}};
    case 'ecog'
        Coords.ECoG= Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}};
        Coords.ECoG_atlas = Location{:,{;'E_atlas_label_Destrieux'}};   
end