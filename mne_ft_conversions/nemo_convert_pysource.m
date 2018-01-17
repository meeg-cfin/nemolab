function source_ft = nemo_convert_pysource(source_mne, fwd_mne)

% NEMO_CONVERT_PYSOURCE converts a source estimate done in MNE-Python and 
% read into MATLAB with mne_read_stc_file to a Fieldtrip source structure.
%
% Input:
% ------
% source_mne : The source structure from MNE-Python as read with 
%              mne_read_stc_file
% fwd_mne : Forward model from MNE-Python as read with 
%           mne_read_forward_solution. Needed for the source space grid
%           positions and transformation matrices.
%
% Output:
% -------
% source_ft : Fieldtrip-like source structure with all relevant information
%             for further usage with plotting functions like ft_sourceplot
%             or nmt_sourceplot.
%
% Author: Britta Westner
%


% initialize
source_ft = [];

% time
source_ft.time = [source_mne.tmin : source_mne.tstep : (source_mne.tmin + source_mne.tstep * (size(source_mne.data,2)-1))];

% positions
source_ft.pos = fwd_mne.src.rr * 1000;  % this is in MEG head space, convert to mm

% inside positions
n_gridpoints = size(source_ft.pos, 1);
inside_idx = source_mne.vertices + 1;  % mne_read_stc_file doesn't seem to correct for Python vs MATLAB indexing
source_ft.inside = zeros(1, n_gridpoints);
source_ft.inside(inside_idx) = 1;
source_ft.inside = logical(source_ft.inside);

% dimensions
% convert to MRI space first b/c axes are aligned here
source_pos_mri = nmt_transform_coord(inv(fwd_mne.mri_head_t.trans), source_ft.pos);
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

end
