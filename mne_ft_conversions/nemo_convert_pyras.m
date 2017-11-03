function pos_nii = nemo_convert_pyras(pos, mri_mgz, mri_nii)

% NEMO_CONVERT_PYRAS converts coordinates in MNE-Python MRI RAS space 
% (Freesurfer RAS coordinates) to nifti MRI voxel space based on the 
% transforms saved in the MRIs.
%
% NOTE: converting the mgz MRI to nifti with Fieldtrip/SPM does not work. 
% Use Freesurfer's mri_convert function for this
% 
% NOTE: Make sure the objects are all in the same units.
%
% Dependencies: NutMEG
%
% Author: Britta Westner

if(~exist('nut_coordtfm'))
    error('Please add NutMEG to your path. Cannot find nut_coordtfm.')
end

% go from Freesurfer RAS space to Freesurfer voxel space:
pos_tmp = nut_coordtfm(pos, inv(mri_mgz.hdr.tkrvox2ras));

% go from Freesurfer voxel space to nifti voxel space:
pos_tmp = nut_coordtfm(pos_tmp, mri_mgz.transform);
pos_nii = nut_coordtfm(pos_tmp, inv(mri_nii.transform));

end
