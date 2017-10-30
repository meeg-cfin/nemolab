% work on getting fieldtrip plotting to work

addpath /home/britta/Documents/git/bw_helper_functions/

figure;
bw_plot3(fwd_model.src.rr, {'b.'}) ; hold on
bw_plot3(fwd_model.source_rr, {'r.', 'markersize', 20});
axis equal

%% let's try to trick Fieldtrip and build a fieldtrip style fwd model

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
gridres = mode(diff(fwd_model.src.rr(:,1)));
dim1 = round((norm(max(lead.pos(1,1)) - min(lead.pos(2,1)))/gridres));
dim2 = round((norm(max(lead.pos(1,2)) - min(lead.pos(2,2)))/gridres));
dim3 = round((norm(max(lead.pos(1,3)) - min(lead.pos(2,3)))/gridres));


%%

ori=svd(lead.pos);
test = lead.pos.*ori(1);
% test = test*ori(2);
% test = test*ori(3);

figure
bw_plot3(test, {'b.'});


%%

% convert to MRI space first b/c axes are aligned here
lp = nut_coordtfm(lead.pos, inv(fwd_model.mri_head_t.trans));

dim(1) = round((max(lp(:,1))-min(lp(:,1))) / gridres);
dim(2) = round((max(lp(:,2))-min(lp(:,2))) / gridres);
dim(3) = round((max(lp(:,3))-min(lp(:,3))) / gridres);

dim = dim+1;


%%
maxp = max(lead.pos); % should be on a plane
gridres = mode(diff(lead.pos(:,1)));

point = zeros(3, 3);
for ii=1:3
    point_idx = dsearchn(lead.pos(:,ii), maxp(ii));
    point(ii,:) = lead.pos(point_idx,:);
end

for ii=1:3
    if ii<3
        dim(ii) = round(norm(point(ii,:)-point(ii+1,:))/gridres);
    else
        dim(ii) = round(norm(point(3,:)-point(1,:))/gridres);
    end
end

dim = sort(dim, 'ascend');
dim = dim+1;

dim(3) = length(lead.pos)/(dim(1)*dim(2))









%%


count=0;
gridres_r = round(gridres, 3);
for ii=1:(length(lead.pos)-1)
    
   step = round(norm(lead.pos(ii,:)-lead.pos(ii+1,:)),3);
   
   if(step ~= gridres_r)
       count=count+1;
       test(count) = ii;       
   end  
    
    
end


%%


lead.dim = [dim1, dim2, dim3]
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
base_path = '/Data2/britta/playing_with_mnedata/';  % in this case, MNE sample data were used

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
% cfg.lcmv.weightnorm    = 'unitnoisegain';
cfg.lcmv.lambda        = '5%';    
source_lcmv            = ft_sourceanalysis(cfg, cov);

%% plot this
source_lcmv.dim = dim;

cfg = [];
cfg.funparameter = 'avg.pow';
ft_sourceplot(cfg, source_lcmv)


%% 





