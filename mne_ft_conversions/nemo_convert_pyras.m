function pos_nii = nemo_convert_pyras(pos, mri_mgz, mri_nii)

% NEMO_CONVERT_PYRAS converts coordinates in MNE-Python MRI RAS space 
% (Freesurfer RAS coordinates) to nifti MRI voxel space based on the 
% transforms saved in the MRIs.
%
% Input:
% ------
% pos : Original positions in RAS surface coordinates
% mri_mgz : mat-file of the original T1.mgz
% mri_nii : mat-file of the T1.mgz converted to nifti (converted with 
%           Freesurfer, see below)
%
% Output:
% -------
% pos_nii : Positions converted to the coordinate frame of the T1.nii
%
% NOTE: converting the mgz MRI to nifti with Fieldtrip/SPM does not work. 
% Use Freesurfer's mri_convert function for this
% 
% NOTE: Make sure the objects are all in the same units.
%
% Dependencies: NutMEG
%
% Author: Britta Westner
%

if(~exist('nmt_transform_coord'))
    error('Please add the folder ./contrib/nutmegtrip of Fieldtrip to your path. Cannot find nmt_transform_coord.')
end

% go from Freesurfer RAS space to Freesurfer voxel space:
pos_tmp = nmt_transform_coord(inv(mri_mgz.hdr.tkrvox2ras), pos);

% go from Freesurfer voxel space to nifti voxel space:
pos_tmp = nmt_transform_coord(mri_mgz.transform, pos_tmp);
pos_nii = nmt_transform_coord(inv(mri_nii.transform), pos_tmp);

end
