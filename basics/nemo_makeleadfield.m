function [grid,vol] = nemo_makeleadfield(cfgnemo)

%% headmodel

switch(cfgnemo.headmodelstrategy)
    case 'singleshell'
        cfg                       = [];
        cfg.grad                  = cfgnemo.grad_mri;
        cfg.feedback              = false;
        cfg.method                = 'singleshell';
        cfg.tissue                = 'brain';
        vol                       = ft_prepare_headmodel(cfg,cfgnemo.bnd(end));
        
        
    case {'dipoli','bemcp'}
        cfg                       = [];
        cfg.grad                  = cfgnemo.grad_mri;
        cfg.feedback              = false;
        cfg.method                = cfgnemo.headmodelstrategy;
        vol                       = ft_prepare_headmodel(cfg,cfgnemo.bnd);

        
    case 'openmeeg'
        cfg                       = [];
        subjId = [cfgnemo.participant '_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer'];
        vol.bnd = cfgnemo.bnd;
        switch(cfgnemo.numlayers)
            case 4
                vol.cond = [0.33 0.0041 1.79 0.33]; % SI units, all 4 layers
                vol.cond = [0.33 0.022 1.79 0.33]; % SI units, all 4 layers % <- from Oostendorp
            case 3
                vol.cond = [0.33 0.0041 0.33]; % SI units, ignore CSF
%                vol.cond = [0.33 0.022 0.33]; % SI units, ignore CSF % <- from Oostendorp
        end
        vol.type = cfgnemo.headmodelstrategy;
        vol.basefile = subjId;
        vol.path = ['./' subjId '/hm/openmeeg_out']; % following files in here can be reused: hm.bin, hm_inv.bin, dsm.bin
        vol = ft_convert_units(vol,'mm');    % Convert bnd to SI units
   case 'simbio'
        cfg                       = [];
        cfg.method = 'simbio';
        subjId = [cfgnemo.participant '_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer'];
        vol.bnd = cfgnemo.bnd;
        switch(cfgnemo.numlayers)
            case 4
                cfg.conductivity = [0.33 0.0041 1.79 0.33]; % SI units, all 4 layers
                cfg.conductivity = [0.33 0.022 1.79 0.33]; % SI units, all 4 layers % <- from Oostendorp
            case 3
                cfg.conductivity = [0.33 0.0041 0.33]; % SI units, ignore CSF
%                vol.cond = [0.33 0.022 0.33]; % SI units, ignore CSF % <- from Oostendorp
        end
        vol.type = cfgnemo.headmodelstrategy;
        vol = ft_convert_units(vol,'mm');    % Convert bnd to SI units
        vol = ft_prepare_headmodel(cfg,vol.bnd);

end

        
%% plotting the headmodel
if(cfgnemo.plotvol)
    for ii=1:length(vol)
        figure
        title(['vol']);
        ft_plot_sens(cfgnemo.grad_mri);
        ft_plot_vol(vol,'facecolor','cortex');
    end
end


%%

% prepare MNI voxel list
cfg                       = [];
cfg.grad                  = cfgnemo.grad_mri; % mm
cfg.vol                   = ft_convert_units(vol,'mm');
cfg.reducerank = 'no';

if(1) % new way: load pre-normalized MRI and determine grid positions from it
    params = load(['s' cfgnemo.participant '_sn.mat']);
    cfg.grid.pos = ft_warp_apply(params, cfgnemo.sourcemodel.pos, 'sn2individual');
else
    % create the subject specific grid, using the template grid that has just been created
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template = cfgnemo.sourcemodel;
    cfg.grid.nonlinear = 'yes'; % use non-linear normalization
    %cfg.grid.resolution = cfgnemo.gridresolution;
    cfg.grid.unit = 'mm',
    cfg.mri            = cfgnemo.mri;
end

if(cfgnemo.VOeyes) % add eyes to "inside" grid
    cfg.grid = ft_prepare_sourcemodel(cfg); % this computes the grid.inside manually
    
    [x,y,z]=meshgrid(-50:50,45:75,-55:-25); % corresponds to eye region in MNI
    eye_xyz = [x(:) y(:) z(:)];
    
    [~,eye_idx]=intersect(cfgnemo.sourcemodel.pos,eye_xyz,'rows');
    cfg.grid.inside(eye_idx) = 1;
end

% cfg.channel = cfgnemo.megchans;

grid               = ft_prepare_leadfield(cfg);

