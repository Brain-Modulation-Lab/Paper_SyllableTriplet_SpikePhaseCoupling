function [Locmap, q_XYZ] = create_LocMap(Location, res_sphere, refMNI, NIFTINAME)
[I1,I2,I3] = ind2sub(refMNI.dim, 1:numel(refMNI.img));
q_XYZ = ea_vox2mm([I1',I2',I3'], refMNI.mat);
[UniqueLocation,nUniqueLocation] = find_uniques(Location);
Idx_Location = rangesearch(UniqueLocation, q_XYZ, res_sphere);
Locmap = cellfun(@(x) sum(nUniqueLocation(x),'omitnan'), Idx_Location);


if ~exist('NIFTINAME','var')
    NIFTINAME = [];
end

if  ~isempty(NIFTINAME)
    refMNI.dt = [16,0];
    refMNI.img = zeros(refMNI.dim);

    % n-map
    LocMap_nii = refMNI;
    LocMap_nii.img = zeros(refMNI.dim);
    LocMap_nii.img(:) = Locmap;
    LocMap_nii.fname = NIFTINAME;

    if exist('ea_space') == 2
        ea_write_nii(LocMap_nii)
    else
        spm_write_vol(LocMap_nii,LocMap_nii.img);
    end
    gzip(LocMap_nii.fname);
    delete(LocMap_nii.fname);
end


end


function [uniques, nuniques] = find_uniques(x)
[uniques,~,ix_uniques] = unique(x,'rows');
nuniques = accumarray(ix_uniques,1)';
end