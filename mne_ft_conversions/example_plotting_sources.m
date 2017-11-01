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
vol_fname = '/sample-5120-bem.fif';
mri_mgz_fname = 'T1.mgz';
mri_nii_fname = 'T1.nii';
fwd_fname = 'sample-fwd.fif';


%% Read the MRIs

mri_mgz = ft_read_mri(fullfile(base_path, mri_mgz_fname));
mri_nii = ft_read_mri(fullfile(base_path, mri_nii_fname));

%% Build a Fieldtrip style foward model

fwd_model = mne_read_forward_solution(fullfile(base_path, fwd_fname), false,  false )
count = 0;
leadfield = cell(1, length(fwd_model.src.rr));

for ii=1:length(fwd_model.src.rr)
    
    if(fwd_model.src.inuse(ii))
        count = count + 1;
        
        leadfield{ii}(:,1) = fwd_model.sol.data(:, count*3-2);   % this is ordered: channels x sourcepos/ori 
        leadfield{ii}(:,2) = fwd_model.sol.data(:, count*3-1);   
        leadfield{ii}(:,3) = fwd_model.sol.data(:, count*3);            
    end     
end

% build the leadfield:
lead.leadfield = leadfield;
lead.pos = fwd_model.src.rr;
lead.inside = logical(fwd_model.src.inuse);
lead.unit = 'm';
lead.label = fwd_model.sol.row_names;
lead.leadfielddimord = '{pos}_chan_ori';

% we need to fake dimensions, too
% convert to MRI space first b/c axes are aligned here
lead_pos_mri = nut_coordtfm(lead.pos, inv(fwd_model.mri_head_t.trans));
gridres = mode(diff(lead_pos_mri(:,1)));

dim(1) = round((max(lead_pos_mri(:,1))-min(lead_pos_mri(:,1))) / gridres);
dim(2) = round((max(lead_pos_mri(:,2))-min(lead_pos_mri(:,2))) / gridres);
dim(3) = round((max(lead_pos_mri(:,3))-min(lead_pos_mri(:,3))) / gridres);

dim = dim+1;
lead.dim = dim;


%% data

% data can be read in with Fieldtrip's preprocessing function:
cfg =[];
cfg.dataset = fullfile(base_path, 'sample_left_vis-epo.fif');
data = ft_preprocessing(cfg);

% select grads only:
cfg = [];
cfg.channel = 'MEGGRAD';
data = ft_selectdata(cfg, data);

%% get data covariance matrix

cfg = [];
cfg.channel          = lead.label;
cfg.covariance       = 'yes';
cfg.covariancewindow = [0.05 0.3];
cov                  = ft_timelockanalysis(cfg, data);  

%% load and convert the volume to head space

% Load the volume:
vol_mne = mne_read_bem_surfaces(fullfile(base_path, vol_fname));

% construct fieldtrip-like structure
vol_ft         = [];
vol_ft.bnd.pos = vol_mne.rr * 1000; % convert to mm
vol_ft.bnd.tri = vol_mne.tris;
vol_ft.type    = 'openmeeg';
vol_ft.cond    = [0.33, 0.0041, 0.33]; % change to whatever was used in MNE

% convert positions from MRI space to MEG head space
vol_ft.bnd.pos = nut_coordtfm(vol_mne.rr, fwd_model.mri_head_t.trans);

%% convert all our inputs to mm

vol_ft = ft_convert_units(vol_ft, 'mm');
lead = ft_convert_units(lead, 'mm');
data.grad = ft_convert_units(data.grad, 'mm');

%% double check: plot everything together

if(1)
    
    figure
    ft_plot_vol(vol_ft, 'facecolor', 'cortex', 'edgecolor', 'None'); alpha 0.2; hold on
    bw_plot3(lead.pos(lead.inside,:), {'r.'});
    bw_plot3(data.grad.chanpos, {'b*'});
    title('Head model, source model, and sensors')
    axis equal
    
end

%% LCMV beamformer

cfg                    = [];
cfg.channel            = cov.label;
cfg.headmodel          = vol_ft;
cfg.method             = 'lcmv';
cfg.grid               = lead;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.reducerank    = 'no' ;
cfg.lcmv.fixedori      = 'yes';
cfg.lcmv.projectnoise  = 'yes';
cfg.lcmv.weightnorm    = 'unitnoisegain';
cfg.lcmv.lambda        = '5%';    
source_lcmv            = ft_sourceanalysis(cfg, cov);

%% plot this without MRI

cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_lcmv)

%% load the fiducials and convert

% load the fiducial file saved from MNE-Python with get_fiducial_info.py
fid_fname = 'test2.mat';
fid_info = load(fullfile(base_path, 'test2'));

% all this is in device coordinates, but we want that thing to go to head coordinates, so:
fid_info_mric = fid_info;

fid_info_mric.fid_coord = nut_coordtfm(fid_info.fid_coord, inv(fwd_model.mri_head_t.trans));
fid_info_mric.fid_coord = fid_info_mric.fid_coord * 1000;
fid_info_mric.fid_coord = nemo_convert_pyras(fid_info_mric.fid_coord, mri_mgz, mri_nii);

%% check whether those fit the MRI

if(1)
    figure;
    ft_plot_ortho(mri_mgz.anatomy, 'transform', mri_mgz.transform, 'style', 'intersect');
    colors = {'b', 'g', 'r'};
    
    hold on
    for ii=1:3
        bw_plot3(fid_info_mric.fid_coord(ii, :), {[colors{ii},'*'], 'markersize', 40}); hold on
        title_maker{ii*2-1} = colors{ii};
        title_maker{ii*2} = fid_info_mric.fid_id{ii};
    end
end
title(sprintf('%s: %s, %s: %s, %s: %s', title_maker{:}))

%% realign the nifti MRI

% Determine the fiducial convention based on the data
switch data.grad.type
    case 'neuromag306'
        head_coord_sys = 'neuromag';
    case {'ctf151', 'ctf275'}
        head_coord_sys = 'ctf';
    otherwise
        error('Head coordinate system type not known, please add.')
end

% Note: this MRI needs to be saved with Freesurfer, Fieldtrip might not work
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = head_coord_sys;
cfg.parameter = 'anatomy';
cfg.viewresult = 'no';
cfg.fiducial.nas = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id, 'nas')), :);
cfg.fiducial.rpa = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id, 'rpa')), :);
cfg.fiducial.lpa = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id, 'lpa')), :);
mri_align = ft_volumerealign(cfg, mri_nii)

%% check whether fiducials fit the aligned MRI

if(1)
    ft_plot_ortho(mri_align.anatomy, 'transform', mri_align.transform, 'style', 'intersect');
    hold on
    for ii=1:3
        bw_plot3(fid_info.fid_coord(ii, :)*1000, {[colors{ii},'*'], 'markersize', 40})
         title_maker{ii*2-1} = colors{ii};
        title_maker{ii*2} = fid_info_mric.fid_id{ii};
    end
    
    title(sprintf('%s: %s, %s: %s, %s: %s', title_maker{:}))

end


%% Interpolate source and plot on the aligned MRI

cfg = [];
cfg.parameter = 'avg.pow';
% cfg.downsample = 2;
source_interp = ft_sourceinterpolate(cfg, source_lcmv, mri_align);

cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_interp);

