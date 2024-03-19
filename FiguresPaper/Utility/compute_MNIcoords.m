function atlas_region = compute_MNIcoords(electrode,cortex_MNI)

% get anatomical information
Vertices = cortex_MNI.Vertices;
if istable(electrode)

    atlas_region = arrayfun(@(x) '', 1:height(electrode), 'uniformoutput', false)';
    for i = find(electrode.type=="ecog")'
        a = project2verts( {electrode{i,{'mni_nonlinear_x', 'mni_nonlinear_y', 'mni_nonlinear_z'}}}, Vertices );
        idx = find(arrayfun(@(x) ismember(a{1}, x.Vertices), cortex_MNI.Atlas(2).Scouts));
        atlas_region{i} = cortex_MNI.Atlas(2).Scouts(idx).Label;
    end
elseif isnumeric(electrode)
    if size(electrode,2) == 3
        atlas_region = cell(1,size(electrode,1));
        for i = 1 : size(electrode,1)
            a = project2verts({electrode(i,:)}, Vertices );
            idx = find(arrayfun(@(x) ismember(a{1}, x.Vertices), cortex_MNI.Atlas(2).Scouts));
            atlas_region{i} = cortex_MNI.Atlas(2).Scouts(idx).Label;
        end
    else
        warning('You need X, Y and Z coordinates.')
    end



end
end