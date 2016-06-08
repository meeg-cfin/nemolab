function grid = nemo_makeleadfield(cfgnemo)
gridresolution = 10

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
        subjId = ['test_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer'];
        vol.bnd = cfgnemo.bnd;
        switch(cfgnemo.numlayers)
            case 4
                vol.cond = [0.33 0.0041 1.79 0.33]; % SI units, all 4 layers
                vol.cond = [0.33 0.022 1.79 0.33]; % SI units, all 4 layers
            case 3
                vol.cond = [0.33 0.0041 0.33]; % SI units, ignore CSF
                vol.cond = [0.33 0.022 0.33]; % SI units, ignore CSF % <- from Oostendorp
        end
        vol.type = cfgnemo.headmodelstrategy;
        vol.basefile = subjId;
        vol.path = ['./' subjId '/hm/openmeeg_out']; % following files in here can be reused: hm.bin, hm_inv.bin, dsm.bin
        vol = ft_convert_units(vol,'mm');    % Convert bnd to SI units

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

switch(cfgnemo.headmodelstrategy)
    case {'singleshell','dipoli','bemcp'}        
        
        % create the subject specific grid, using the template grid that has just been created
        cfg.grid.warpmni   = 'yes';
        load standard_sourcemodel3d10mm; % loads in sourcemodel (i.e., MNI voxel grid)
        sourcemodel  = ft_convert_units(sourcemodel,'mm');
        cfg.grid.template = sourcemodel;
        cfg.grid.nonlinear = 'yes'; % use non-linear normalization
        cfg.mri            = cfgnemo.mri;
        clear mri
        cfg.unit = 'mm',
        
        cfg.vol                   = ft_convert_units(vol,'mm');
        cfg.megchans = cfgnemo.megchans;
        cfg.reducerank = 'no';
        cfg.grid.resolution = gridresolution;
        grid               = ft_prepare_leadfield(cfg);
        

    case 'openmeeg'
        switch(cfgnemo.voxelgridtype)
            case 'subjectspace'
                cfg.vol                   = ft_convert_units(vol,'mm');
                cfg.grid.resolution       = 10; % use 10mm for quick look; 5 for final computation
                cfg.grid.unit             = 'mm';
                cfg.megchans = megchans;
                cfg.reducerank = 'no';
                grid=ft_prepare_leadfield(cfg);
            case 'mni-ft'
                % create the subject specific grid, using the template grid that has just been created
                cfg.grid.warpmni   = 'yes';
                %standard_sourcemodel3d5mm;  % in cm!!!!
                load standard_sourcemodel3d10mm; % loads in sourcemodel (i.e., MNI voxel grid)
                %            load template_grid; sourcemodel = template_grid;
                sourcemodel  = ft_convert_units(sourcemodel,'mm');
                cfg.grid.template = sourcemodel;
                cfg.grid.nonlinear = 'yes'; % use non-linear normalization
                cfg.mri            = cfgnemo.mri;
                cfg.grid.unit = 'mm',
                %    gridmni               = ft_prepare_sourcemodel(cfg);
                
                cfg.vol                   = ft_convert_units(vol,'mm');
                cfg.megchans = cfgnemo.megchans;
                cfg.reducerank = 'no';
                cfg.grid.resolution = gridresolution;
                cfg.grad                  = ft_convert_units(cfgnemo.grad_mri,'mm'); % SI units
                cfg.megchans = cfgnemo.megchans;
                cfg.reducerank = 'no';
                grid               = ft_prepare_leadfield(cfg);
                
                
                
            case 'mni-nm'
                % TEST NEEDED: does this work at all? is it compatible with NMT plotting?
                load template_grid
                grid = template_grid;
                grid.pos = nut_mni2mri(template_grid.pos);
        end
end

cfg = [];