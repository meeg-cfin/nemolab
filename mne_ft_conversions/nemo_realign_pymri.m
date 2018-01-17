function mri_aligned = nemo_realign_pymri(mri_nii, fid_info, sens_type, ras2meg_tfm, mri_mgz)

% NEMO_REALIGN_PYMRI realigns a nifti MRI to MEG head space as defined by
% fidcuials imported from MNE-Python.
%
% Input:
% ------
% mri_nii : nifti version of the original T1.mgz as mat-file, converted 
%           with Freesurfer
% fid_info : Fidcuial information containing the fiducials positions in 
%            MEG head space. Saved from Python with the function
%            get_fiducial_info.py
% sens_type : information about the sensor type, e.g. 'neuromag'. For more
%             information see ft_volumerealign : cfg.coordsys option.
% ras2meg_tfm : Transform between RAS surface coordinates and MEG head 
%               space, e.g. from forwardmodel.mri_head_t.trans
% mri_mgz : mat-file of the original T1.mgz
%
%
% Output:
% -------
% mri_aligned : mat-file containing the MRI aligned to MEG head space
%
% Dependencies: Fieldtrip, including NutMEGtrip
%
% Author: Britta Westner
%


if(~exist('ft_volumerealign'))
    error('Please add FieldTrip to your path. Cannot find ft_volumerealign.')
end

% convert fiducials from MEG head space to nifti space:
fid_info.fid_coord = nmt_transform_coord(inv(ras2meg_tfm), fid_info.fid_coord);
fid_info.fid_coord = fid_info.fid_coord * 1000;
fid_info.fid_coord = nemo_convert_pyras(fid_info.fid_coord, mri_mgz, mri_nii);

% realign the MRI to MEG head space based on the fiducials:
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = sens_type;
cfg.parameter = 'anatomy';
cfg.viewresult = 'no';
cfg.fiducial.nas = fid_info.fid_coord(find(strcmp(fid_info.fid_id, 'nas')), :);
cfg.fiducial.rpa = fid_info.fid_coord(find(strcmp(fid_info.fid_id, 'rpa')), :);
cfg.fiducial.lpa = fid_info.fid_coord(find(strcmp(fid_info.fid_id, 'lpa')), :);

mri_aligned = ft_volumerealign(cfg, mri_nii);

end
