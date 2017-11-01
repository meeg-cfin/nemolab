% work on getting fieldtrip plotting to work

addpath /home/britta/Documents/git/bw_helper_functions/

% DATA
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



%% let's try to trick Fieldtrip and build a fieldtrip style fwd model
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
datapath = '/Data2/britta/playing_with_mnedata/sample_left_vis-epo.fif'
cfg =[];
cfg.dataset = datapath;
data = ft_preprocessing(cfg);

% select grads only:
cfg = [];
cfg.channel = 'MEGGRAD';
data = ft_selectdata(cfg, data);

%% use this to source reconstruct

% create Fieldtrip covariance matrix from MNE data
cfg = [];
cfg.channel          = lead.label;
cfg.covariance       = 'yes';
cfg.covariancewindow = [0.05 0.3];
cov                  = ft_timelockanalysis(cfg, data);   % dummy

%% load and convert the volume to head space

% Filenames, should be relative to base_path
vol_fname = '/sample-5120-bem.fif';
vol_mne = mne_read_bem_surfaces(fullfile(base_path, vol_fname));

% construct fieldtrip structure
vol_ft         = [];
vol_ft.bnd.pos = vol_mne.rr * 1000; % convert to mm
vol_ft.bnd.tri = vol_mne.tris;
vol_ft.type    = 'openmeeg';
vol_ft.cond    = [0.33, 0.0041, 0.33]; % change to whatever was used in MNE

% replace positions with pos in MEG head space
vol_ft.bnd.pos = nut_coordtfm(vol_mne.rr, fwd_model.mri_head_t.trans);

%% convert all to mm

vol_ft = ft_convert_units(vol_ft, 'mm');
lead = ft_convert_units(lead, 'mm');
data.grad = ft_convert_units(data.grad, 'mm');
%% double check
if(1)
    
    figure
    ft_plot_vol(vol_ft, 'facecolor', 'cortex', 'edgecolor', 'None'); alpha 0.2; hold on
    bw_plot3(lead.pos(lead.inside,:), {'r.'});
    bw_plot3(data.grad.chanpos, {'b*'});
    title('Head model and forward model')
    axis equal
    
end

%% convert all to nifti space
ras2meg_mm = fwd_model.mri_head_t.trans;
ras2meg_mm(1:3, 4) = ras2meg_mm(1:3, 4) * 1000;
 
vol_mri = vol_ft;
vol_mri.bnd.pos = nemo_convert_pyras(vol_mne.rr * 1000, mri_mgz, mri_nii);

lead_mri = lead;
lead_mri.pos = nut_coordtfm(lead_mri.pos, inv(ras2meg_mm));
lead_mri.pos = nemo_convert_pyras(lead_mri.pos, mri_mgz, mri_nii);

data_mri = data;
data_mri.grad.chanpos = nut_coordtfm(data_mri.grad.chanpos, inv(ras2meg_mm));
data_mri.grad.chanpos = nemo_convert_pyras(data_mri.grad.chanpos, mri_mgz, mri_nii);
data_mri.grad.coilpos = nut_coordtfm(data_mri.grad.coilpos, inv(ras2meg_mm));
data_mri.grad.coilpos = nemo_convert_pyras(data_mri.grad.coilpos, mri_mgz, mri_nii);

%% double check II

if(1)
    
    figure
    ft_plot_vol(vol_mri, 'facecolor', 'cortex', 'edgecolor', 'None'); alpha 0.2; hold on
    bw_plot3(lead_mri.pos(lead_mri.inside,:), {'r.'});
    bw_plot3(data_mri.grad.chanpos, {'b*'});
    title('Head model and forward model')
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

%% plot this

cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_lcmv)


%% nutmegtrip

fieldtrippathnmt = '/Data/MATLAB/dev/fieldtrip'%'/home/britta/Documents/Matlab_toolboxes/fieldtrip-20170606';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));

cfg=[];
% cfg.zlim = [-3 3];
% inside_idx = find(source_tf.inside);
cfg.mripath = fullfile(base_path, mri_nii_fname);
cfg.funparameter = 'avg.pow';
% cfg.evokedoverlay = 1;
% cfg.plottype = 'tf';  % ts if you want to view as time series instead
cfg.atlas = ft_read_atlas(fullfile(fieldtrippathnmt, '/template/atlas/aal/ROI_MNI_V4.nii')); % AAL atlas
nmt_sourceplot(cfg,ft_convert_units(source_lcmv,'mm'));

%% plot this on the mgz mri
mri_mgz2 = mri_mgz;
mri_mgz2.transform = mri_mgz.transform*inv(mri_mgz.hdr.tkrvox2ras)*(fwd_model.mri_head_t.trans);

cfg = [];
cfg.parameter = 'avg.pow';
cfg.downsample = 5;
source_interp = ft_sourceinterpolate(cfg, source_lcmv, mri_mgz);

%% 

% source_pos = nut_coordtfm(source_interp.pos, inv(fwd_model.mri_head_t.trans));
% source_pos = nut_coordtfm(source_interp.pos, inv(mri_mgz.transform));
% source_cp = source_interp;
% source_cp.pos = source_pos;
% source_cp.transform = inv(mri_mgz.transform);

cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_interp);

%% try plotting with nutmegtrip instead


% convert to nii
% source_pos = nut_coordtfm(source_lcmv.pos, inv(fwd_model.mri_head_t.trans));
% source_pos = nemo_convert_pyras(source_pos, mri_mgz, mri_nii);

source_cp = source_lcmv;
% source_cp.pos = source_pos;




fieldtrippathnmt = '/Data/MATLAB/dev/fieldtrip'%'/home/britta/Documents/Matlab_toolboxes/fieldtrip-20170606';
addpath(fullfile(fieldtrippathnmt, 'contrib/nutmegtrip'));

fname = 'mri_tfm.nii';
nii_tfm = inv(meg_head_space_tfm)*mri_mgz.hdr.tkrvox2ras * inv(mri_mgz.transform) * mri_nii.transform;
mri_nii_new = mri_nii;
mri_nii_new.transform = nii_tfm;
ft_write_mri(fullfile(base_path, fname), mri_nii_new, 'dataformat', 'nifti');


mri_path = fullfile(base_path, fname);


cfg=[];
cfg.zlim = [-3 3];
% inside_idx = find(source_tf.inside);
cfg.mripath = mri_path;
cfg.funparameter = 'avg.pow';
% cfg.evokedoverlay = 1;
% cfg.plottype = 'tf';  % ts if you want to view as time series instead
cfg.atlas = ft_read_atlas(fullfile(fieldtrippathnmt, '/template/atlas/aal/ROI_MNI_V4.nii')); % AAL atlas
nmt_sourceplot(cfg,ft_convert_units(source_cp,'mm'));

%% convert mri to head space

% go from MRI voxel to MEG head space:
% 0) apply mri transform
% 1) apply tkrvox2ras
% 2) apply inv(mri_head_t.trans)

meg_head_space_tfm = fwd_model.mri_head_t.trans;
meg_head_space_tfm(1:3, 4) = meg_head_space_tfm(1:3, 4)*1000;
source_lcmv.transform = (meg_head_space_tfm);

tfm = ( (meg_head_space_tfm)*mri_mgz.hdr.tkrvox2ras) %* mri_mgz.transform;

mri_mgz_meg = mri_mgz;
mri_mgz_meg.transform = tfm;
%%
cfg = [];
cfg.parameter = 'avg.pow';
% cfg.downsample = 2;
source_interp = ft_sourceinterpolate(cfg, source_lcmv, mri_align);


cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_interp);

%% load fids and convert

fid_info = load(fullfile(base_path, 'test2'));

% all this is in device coordinates, but we want that thing to go to head 
% coordinates, so:
fid_info_mric = fid_info;

fid_info_mric.fid_coord = nut_coordtfm(fid_info.fid_coord, inv(fwd_model.mri_head_t.trans));
fid_info_mric.fid_coord = fid_info_mric.fid_coord * 1000;
fid_info_mric.fid_coord = nemo_convert_pyras(fid_info_mric.fid_coord, mri_mgz, mri_nii);


%% check whether those fit the MRI

ft_plot_ortho(mri_mgz.anatomy, 'transform', mri_mgz.transform, 'style', 'intersect');
% ft_sourceplot([], mri_mgz);
hold on
for ii=1:3
    bw_plot3(fid_info_mric.fid_coord(ii, :), {[colors{ii},'*'], 'markersize', 40}); hold on
end

%% and double check:

figure; hold all
bw_plot3(data.grad.chanpos);
colors = {'b', 'r', 'g'}
for ii=1:3
    
    bw_plot3(fid_info_mric.fid_coord(ii, :), {[colors{ii},'*'], 'markersize', 20})
    title_maker{ii*2-1} = colors{ii};
    title_maker{ii*2} = fid_info_mric.fid_id{ii};
 
end
title(sprintf('%s: %s, %s: %s, %s: %s', title_maker{:}))

%% realign nifti

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'neuromag';
cfg.parameter = 'anatomy';
cfg.viewresult = 'no';
cfg.fiducial.nas = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id,'nas')), :);
cfg.fiducial.rpa = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id,'rpa')), :);
cfg.fiducial.lpa = fid_info_mric.fid_coord(find(strcmp(fid_info_mric.fid_id,'lpa')), :);
mri_align = ft_volumerealign(cfg, mri_nii)

%% fiducials fit the nifti
ft_plot_ortho(mri_align.anatomy, 'transform', mri_align.transform, 'style', 'intersect');

% ft_sourceplot([], mri_align);
hold on
for ii=1:3
    bw_plot3(fid_info.fid_coord(ii, :)*1000, {[colors{ii},'*'], 'markersize', 40})
end
