%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Converting MNE-Python source data to plot them with NutMEGtrip      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script loads the output of the MNE-Python lcmv() and plots it with
% Nutmegtrip
%
% Author: Britta Westner


% Path to MNE-Python output folder
base_path = '/path/to/data';  % in this case, MNE sample data were used

% source estimate filename from MNE-Python, saved with source_est.save()
% NOTE: this script assumes the source reconstruction of evoked data, i.e. one time
% series per voxel
source_fname = 'source_est-vl.stc';

% Files needed for transformations:
mri_mgz_fname = 'T1.mgz';  % .mgz MRI used in Freesurfer
mri_nii_fname = 'T1.nii';  % .nii version of the above MRI, transformed with Freesurfer
fwd_fname = 'sample_2-fwd.fif';  % forward model for transform matrices and source grid


%% read the source space estimate and the forward model

source_mne = mne_read_stc_file(fullfile(base_path, source_fname));
fwd_model = mne_read_forward_solution(fullfile(base_path, fwd_fname), false,  false )


%% make the source estimate a Fieldtrip structure

source_ft = [];

% time
source_ft.time = [source_mne.tmin : source_mne.tstep : (source_mne.tmin + source_mne.tstep * (size(source_mne.data,2)-1))];

% positions 
source_ft.pos = fwd_model.src.rr * 1000;  % this is in MEG head space, convert to mm

% inside positions
n_gridpoints = size(source_ft.pos, 1);
inside_idx = source_mne.vertices + 1;  % mne_read_stc_file doesn't seem to correct for Python vs MATLAB indexing
source_ft.inside = zeros(1, n_gridpoints);
source_ft.inside(inside_idx) = 1;
source_ft.inside = logical(source_ft.inside);

% dimensions
% convert to MRI space first b/c axes are aligned here
source_pos_mri = nut_coordtfm(source_ft.pos, inv(fwd_model.mri_head_t.trans));
gridres = mode(diff(source_pos_mri(:,1)));

dim(1) = round((max(source_pos_mri(:,1))-min(source_pos_mri(:,1))) / gridres);
dim(2) = round((max(source_pos_mri(:,2))-min(source_pos_mri(:,2))) / gridres);
dim(3) = round((max(source_pos_mri(:,3))-min(source_pos_mri(:,3))) / gridres);
dim = dim+1;
source_ft.dim = dim;

% bookkeeping etc.
source_ft.method = 'average';
source_ft.cfg.toolbox = 'MNE-Python import';
source_ft.cfg.gridres = gridres;  % might be interesting to have
source_ft.cfg.units = 'mm';

% dimord
source_ft.avg.filterdimord = '{pos}_ori_chan';  % doesn't make too much sense but that's what Fieldtrip inserts

% moments:
source_ft.avg.mom = cell(n_gridpoints, 1);

for ii=1:length(inside_idx)    
    source_ft.avg.mom{inside_idx(ii)} = source_mne.data(ii, :);
end

% calculate power estimate
source_ft.avg.pow = NaN(n_gridpoints, 1);

for ii=1:length(inside_idx)    
   source_ft.avg.pow(inside_idx(ii)) = sum(source_mne.data(ii, :).^2); 
end
    
%% Read the MRIs

mri_mgz = ft_read_mri(fullfile(base_path, mri_mgz_fname));
mri_nii = ft_read_mri(fullfile(base_path, mri_nii_fname));

%% Convert the positions to nifti

% Positions are in RAS mgz MRI space and need to go to nifti space
ras2meg = fwd_model.mri_head_t.trans;
ras2meg(1:3, 4) = ras2meg(1:3, 4) * 1000;  % convert to mm
source_pos = nut_coordtfm(source_ft.pos, inv(ras2meg));
source_pos = nemo_convert_pyras(source_pos, mri_mgz, mri_nii);
% and go back to common space due to nii transform (not needed with FT plotting)
source_pos = nut_coordtfm(source_pos, mri_nii.transform);

% make a copy of the source to prevent failures with multiple runs of the 
% same script and insert converted positions
source_nii = source_ft;
source_nii.pos = round(source_pos, 3); % rounding prevents mode() failure

%% plot power with NutMEGtrip

% path to Nutmeg
fieldtrippathnmt = '/Data/MATLAB/dev/fieldtrip';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));

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
