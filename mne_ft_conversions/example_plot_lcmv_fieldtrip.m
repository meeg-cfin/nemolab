%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Converting MNE-Python source data to plot them with Fieldtrip      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script loads the output of the MNE-Python lcmv() and plots it with
% Fieldtrip (in MEG head space)
%
% Author: Britta Westner

%% specifications and path names

% TOOLBOXES
% Fieldtrip
fieldtrip_path = '/path/to/fieldtrip';

% DATA
% Path to MNE-Python output folder
base_path = '/path/to/data';

% source estimate filename from MNE-Python, saved with source_est.save()
% NOTE: this script assumes the source reconstruction of evoked data, i.e. one time
% series per voxel
source_fname = 'source_est-vl.stc';

% Files needed for transformations:
mri_mgz_fname = 'T1.mgz';  % .mgz MRI used in Freesurfer
mri_nii_fname = 'T1.nii';  % .nii version of the above MRI, transformed with Freesurfer
fwd_fname = 'sample_2-fwd.fif';  % forward model for transform matrices and source grid
fid_fname = 'test2.mat';  % filename with fiducial information. saved from Python with get_fiducial_info.py

% either one of the following:
data_fname = 'sample_left_vis-epo.fif';  % needed to get data type. 
% if the sensor type is known, it can also be specified directly here:
if(1); sens_type = 'neuromag'; end  % !! sens_type overrules data_fname

%% adding toolboxes to path
% handling Fieldtrip paths and preventing multiple paths from being added:

try 
    ft_defaults    
catch
    warning('Fieldtrip is not on your path yet, adding it.');
    addpath(fieldtrip_path)
    ft_defaults
end

[ft_ver, ft_path] = ft_version;
display(sprintf('You are using Fieldtrip on path %s', ft_path));

% path to Nutmegtrip and MNE externals - Fieldtrip itself not needed here
fieldtrippathnmt = '/path/to/fieldtrip';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));
addpath(fullfile(fieldtrippathnmt, 'external/mne'));

%% read the source space estimate and the forward model

source_mne = mne_read_stc_file(fullfile(base_path, source_fname));
fwd_model = mne_read_forward_solution(fullfile(base_path, fwd_fname), false,  false )

% get the transform out of the fwd_model:
ras2meg_tfm = fwd_model.mri_head_t.trans;

%% make the source estimate a Fieldtrip structure

source_ft = nemo_convert_pysource(source_mne, fwd_model);

%% Read the MRIs

mri_mgz = ft_read_mri(fullfile(base_path, mri_mgz_fname));
mri_nii = ft_read_mri(fullfile(base_path, mri_nii_fname));

%% realign the nifti MRI to head space

% load the fiducial file saved from MNE-Python with get_fiducial_info.py
fid_info = load(fullfile(base_path, fid_fname));

% Determine the fiducial convention based on the data or string input
if(exist('sens_type', 'var'))
    
    head_coord_sys = sens_type;
    
elseif(exist('data_fname', 'var'))    
    
    cfg =[];
    cfg.dataset = fullfile(base_path, data_fname);
    data = ft_preprocessing(cfg);

    switch data.grad.type
        case 'neuromag306'
            head_coord_sys = 'neuromag';
        case {'ctf151', 'ctf275'}
            head_coord_sys = 'ctf';
        otherwise
            error('Head coordinate system type not known, please add.')
    end
else
    error('Need something to infer sensor type from, either data structure or sens_type specification.')
end

% realign the nifti:
mri_aligned = nemo_realign_pymri(mri_nii, fid_info, head_coord_sys, ras2meg_tfm, mri_mgz);

%% Interpolate source and plot on the aligned MRI

cfg = [];
cfg.parameter = 'avg.pow';
% cfg.downsample = 3;  % for faster performance during tests
source_interp = ft_sourceinterpolate(cfg, source_ft, mri_aligned);

cfg = [];
cfg.funparameter = 'avg.pow';
cfg.method = 'slice';
ft_sourceplot(cfg, source_interp);
