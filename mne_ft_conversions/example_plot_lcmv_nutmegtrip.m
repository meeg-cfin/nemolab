%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Converting MNE-Python source data to plot them with NutMEGtrip      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script loads the output of the MNE-Python lcmv() and plots it with
% Nutmegtrip (in nifti space)
%
% Author: Britta Westner

%% specifications and path names

% TOOLBOX
% path to Nutmegtrip and MNE externals - Fieldtrip itself not needed here
fieldtrippathnmt = '/path/to/fieldtrip';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));
addpath(fullfile(fieldtrippathnmt, 'external/mne'));

% DATA
% Path to MNE-Python output folder
base_path = '/path/to/data';  % in this case, MNE sample data were used

% source estimate filename from MNE-Python, saved with source_est.save()
% NOTE: this script assumes the source reconstruction of evoked data, i.e. one time
% series per voxel
source_fname = 'source_est-vl.stc';

% Files needed for transformations:
mri_mgz_fname = 'T1.mgz';  % .mgz MRI used in Freesurfer
mri_nii_fname = 'T1.nii';  % .nii version of the above MRI, transformed with Freesurfer (mri_convert)
fwd_fname = 'sample_2-fwd.fif';  % forward model for transform matrices and source grid

%% read the source space estimate and the forward model

source_mne = mne_read_stc_file(fullfile(base_path, source_fname));
fwd_model = mne_read_forward_solution(fullfile(base_path, fwd_fname), false,  false )

%% make the source estimate a Fieldtrip structure

source_ft = nemo_convert_pysource(source_mne, fwd_model);

%% Read the MRIs

mri_mgz = ft_read_mri(fullfile(base_path, mri_mgz_fname));
mri_nii = ft_read_mri(fullfile(base_path, mri_nii_fname));

%% Convert the positions to nifti

% Positions are in RAS mgz MRI space and need to go to nifti space
ras2meg = fwd_model.mri_head_t.trans;
ras2meg(1:3, 4) = ras2meg(1:3, 4) * 1000;  % convert to mm
source_pos = nmt_transform_coord(inv(ras2meg), source_ft.pos);
source_pos = nemo_convert_pyras(source_pos, mri_mgz, mri_nii);
% and go back to common space due to nii transform (not needed with FT plotting)
source_pos = nmt_transform_coord(mri_nii.transform, source_pos);

% make a copy of the source to prevent failures with multiple runs of the 
% same script and insert converted positions
source_nii = source_ft;
source_nii.pos = round(source_pos, 3); % rounding prevents mode() failure

%% plot power with NutMEGtrip

% plot the power estimate
cfg=[];
cfg.mripath = mri_nii_fname;
cfg.funparameter = 'avg.pow';
nmt_sourceplot(cfg, source_nii);

%% plot time courses with NutMEGtrip

% plot the time courses
cfg=[];
cfg.mripath = mri_nii_fname;
cfg.funparameter = 'avg.mom';
nmt_sourceplot(cfg, source_nii);
