%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Converting MNE-Python data to plot with Fieldtrip             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the reads MNE-Python input (head model and forward model)
% and runs an LCMV beamformer on it to finally plot it with Fieldtrip functions.
%
% Author: Britta Westner


% Path to MNE-Python subject folder
base_path = '/Data2/britta/playing_with_mnedata';  % in this case, MNE sample data were used

% Filenames, should be relative to base_path
% vol_fname = '/sample-5120-bem.fif';
mri_mgz_fname = 'T1.mgz';
mri_nii_fname = 'T1.nii';
fwd_fname = 'sample_2-fwd.fif';
source_fname = 'source_est-vl.stc';

%% read the source space estimate and the forward model

source_mne = mne_read_stc_file(source_fname);
fwd_model = mne_read_forward_solution(fullfile(base_path, fwd_fname), false,  false )


%% make this a Fieldtrip structure

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

ras2meg = fwd_model.mri_head_t.trans;
ras2meg(1:3, 4) = ras2meg(1:3, 4) * 1000;  % convert to mm
source_pos = nut_coordtfm(source_ft.pos, inv(ras2meg));
source_pos = nemo_convert_pyras(source_pos, mri_mgz, mri_nii);
% and go back to common space due to nii transform (not needed with FT plotting)
source_pos = nut_coordtfm(source_pos, mri_nii.transform);

source_ft.pos = round(source_pos, 3);

%% plot power with NutMEGtrip

fieldtrippathnmt = '/Data/MATLAB/dev/fieldtrip';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));

cfg=[];
cfg.mripath = mri_nii_fname;
cfg.funparameter = 'avg.pow';
nmt_sourceplot(cfg,ft_convert_units(source_ft,'mm'));

%% plot time courses with NutMEGtrip

cfg=[];
cfg.mripath = mri_nii_fname;
cfg.funparameter = 'avg.mom';
nmt_sourceplot(cfg,ft_convert_units(source_ft,'mm'));
