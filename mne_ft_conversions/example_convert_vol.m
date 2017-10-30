%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Converting MNE-Python meshes to nifti coordinate frame           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the NEMO_CONVERT_PYRAS function to get MNE-Python BEM
% meshes to the nifti coordinate frame
%
% Author: Britta Westner

%% specification of paths and filenames

% TOOLBOXES
% Fieldtrip
fieldtrip_path = '/path/to/fieldtrip';
% NutMEG
nutmeg_path = '/path/to/nutmeg';

% DATA
% Path to MNE-Python subject folder
base_path = '/path/to/data/';  % in this case, MNE sample data were used

% Filenames, should be relative to base_path
vol_fname = '/sample-5120-bem.fif';
mri_mgz_fname = 'T1.mgz';
mri_nii_fname = 'T1.nii';
fwd_fname = 'sample-fwd.fif';

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

% handling NutMEG:
addpath(nutmeg_path)

%% Loading the BEM and constructing a Fieldtrip-like structure from it

vol_mne = mne_read_bem_surfaces(fullfile(base_path, vol_fname));

% construct fieldtrip structure
vol_ft         = [];
vol_ft.bnd.pos = vol_mne.rr * 1000; % convert to mm
vol_ft.bnd.tri = vol_mne.tris;
vol_ft.type    = 'openmeeg';
vol_ft.cond    = [0.33, 0.0041, 0.33]; % change to whatever was 

%% Read the MRIs

mri_mgz = ft_read_mri(fullfile(base_path, mri_mgz_fname));
mri_nii = ft_read_mri(fullfile(base_path, mri_nii_fname));

%% convert the meshes to nifti

% input the original positions, but converted to mm
pos_tmp = nemo_convert_pyras(vol_mne.rr * 1000, mri_mgz, mri_nii);

vol_nii = vol_ft;
vol_nii.bnd.pos = pos_tmp;

% plot the new mesh together with the nifti
figure;
ft_plot_ortho(mri_nii.anatomy, 'style', 'intersect'); hold on
ft_plot_mesh(vol_nii.bnd, 'facecolor', 'magenta'); 

