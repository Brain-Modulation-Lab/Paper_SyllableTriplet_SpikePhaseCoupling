function NodeMap = LocMap2NodeMap(LocMapPATH,flag_reduced)

if isstring(LocMapPATH) || ischar(LocMapPATH)
    if ~contains(LocMapPATH,'.nii.gz')
        LocMapPATH = strrep(LocMapPATH,'.nii','.nii.gz');
    end
    hdr = spm_vol(LocMapPATH);
    [img, xyz] = spm_read_vols(hdr);
    values = img(:);
else
    warning('error: put a filepath (.nii file)')
end

if ~exist('flag_reduced','var')
    flag_reduced = false;
end


NodeMap = table();
NodeMap.X = xyz(1,:)';
NodeMap.Y = xyz(2,:)';
NodeMap.Z = xyz(3,:)';
NodeMap.Color = values;
NodeMap.Size =  values*3;
NodeMap.Color(isnan(NodeMap.Color)) = 0;
NodeMap.Size(isnan(NodeMap.Size)) = 0;

if ~flag_reduced
    writetable(NodeMap,strrep(LocMapPATH,'.nii.gz','.txt'))
    moveFile(strrep(LocMapPATH,'.nii.gz','.txt'), strrep(LocMapPATH,'.nii.gz','.node'));
    else
    NodeMap(NodeMap.Color == 0,:) = [];
    writetable(NodeMap,strrep(LocMapPATH,'.nii.gz','_reduced.txt'))
    moveFile(strrep(LocMapPATH,'.nii.gz','_reduced.txt'), strrep(LocMapPATH,'.nii.gz','_reduced.node'));
end


