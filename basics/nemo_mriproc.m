
plotvol = true;

%% reslice mri (if not already done)
% reslice_nii('SSD.nii','sSSD.nii',1); % reslice into isotropic 1 mm voxels
%
%% COREGISTER in NUTMEG
% nm
%%% then:
% global nuts; coreg = nuts.coreg; save sSSDcoreg coreg
load(['s' cfgnemo.participant 'coreg.mat']);

switch(data.grad.type)
    case {'bti148','bti248'}
        origcoordsys = 'bti';
    case 'neuromag306'
        origcoordsys = 'neuromag';
    otherwise
        error('Which fiducial convention do you have? Add it to the cases above!')
end
[coreg.mri2meg_tfm,coordsys] = ft_headcoordinates(coreg.fiducials_mri_mm(3,:),coreg.fiducials_mri_mm(1,:),coreg.fiducials_mri_mm(2,:),origcoordsys);
coreg.meg2mri_tfm = inv(coreg.mri2meg_tfm);

%% transform MEG to MRI space
grad = data.grad;  % NOTE data.hdr.grad does not contain correct information if synthetic gradient has been manipulated above
grad_mm = ft_convert_units(grad,'mm'); % transforms the grad units (m) to the same than mri (mm)
grad_mri = ft_transform_sens(coreg.meg2mri_tfm, grad_mm); % transforms the grad coordinates to mri coordinates
grad_mri.coordsys = 'spm';


%% load mri
basepath = pwd;
mripath = [basepath '/s' cfgnemo.participant '.nii'];
normmripath = [basepath '/ws' cfgnemo.participant '.nii'];
mri = ft_read_mri(mripath,'dataformat','nifti');
mri.coordsys = 'spm';


%% segmentation (runs only if seg_*.mat is not already present)
if(exist([basepath '/seg_' cfgnemo.segmethod '.mat'],'file'))
    load([basepath '/seg_' cfgnemo.segmethod '.mat']);
else
    switch cfgnemo.segmethod
        case 'spm8newseg'
            seg = create_bem_segmentation('t1filename',mripath);
            seg.transform = mri.transform;
            seg.dim = mri.dim;
            seg.unit = mri.unit;
        case 'ftvolseg'
            cfg = [];
            mri.coordsys = 'ras';
            cfg.output = {'gray','white','csf','skull','scalp'};
            seg = ft_volumesegment(cfg,mri);
        otherwise
            error('unknown segmentation method requested')
    end
    save([basepath '/seg_' cfgnemo.segmethod '.mat'],'seg');
end

%%
if(exist([basepath '/bnd_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer.mat'],'file'))
    load([basepath '/bnd_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer.mat']);
else
    cfg = [];
    switch(cfgnemo.numlayers)
        case 4
            cfg.tissue = {'scalp', 'skull', 'csf', 'brain'};
        case 3
            cfg.tissue = {'scalp', 'skull', 'brain'};
    end
    
    switch(cfgnemo.headmodelstrategy)
        case 'simbio'
            cfg.method = 'hexahedral';
            cfg.numvertices = [200 400 800 800];
            bnd = ft_prepare_mesh(cfg, seg);
        otherwise
            cfg.method = 'iso2mesh';    % 'projectmesh';
            cfg.numvertices = 10000;    % We'll decimate later
            bnd = ft_prepare_mesh(cfg, seg);
            % Mesh repairs - Not yet implemented in FT
            %    targetsize = [500 1000 1500 1500]; % final mesh size desired for layers (in order given above)
            targetsize = [1000 1000 1000 1000]; % final mesh size desired for layers (in order given above)
            
            % decimate, check, and repair individual meshes using iso2mesh
            for ii = 1:length(bnd)
                [bnd(ii).pos, bnd(ii).tri] = meshresample(bnd(ii).pos, bnd(ii).tri, targetsize(ii)/size(bnd(ii).pos,1));
                [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'dup');
                [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'isolated');
                [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'deep');
                [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'meshfix');
            end
            % Ensure no overlaps
            bnd = nemo_decouplesurf(bnd);    % decouplesurf is an unimplemented subfunction temporarily stashed in prepare_mesh_segmentation
            
            if(cfgnemo.plotvol)
                figure;
                ft_plot_sens(grad_mri);
                ft_plot_mesh(bnd(1), 'facealpha', 0.25), hold on
                ft_plot_mesh(bnd(2), 'facealpha', 0.25, 'facecolor', 'blue')
                ft_plot_mesh(bnd(3), 'facealpha', 0.25, 'facecolor', 'yellow')
                if(length(bnd)>=4)
                    ft_plot_mesh(bnd(4), 'facecolor', 'red')
                end
                ft_plot_ortho(mri.anatomy,'transform',mri.transform,'style','intersect')
                % plot3(grad_mri.chanpos(:,1),grad_mri.chanpos(:,2),grad_mri.chanpos(:,3),'go')
            end
    end
    
    clear seg
    save([basepath '/bnd_' cfgnemo.segmethod '_' num2str(cfgnemo.numlayers) 'layer.mat'],'bnd');
    
    
    
end

cfgnemo.mri = mri;

%% Spatial normalization via SPM
%
% if using SPM12, SPM12/toolbox/OldNorm needed in path
% TODO: try out SPM12's DARTEL normalization

if(~exist(normmripath,'file'))
    cfg = [];
    tmpcfg.nonlinear = 'yes';
    spmdir = fileparts(which ('spm'));
    tmpcfg.spmversion = spm('ver');
    switch(spm('ver'))
        case 'SPM8'
            tmpcfg.template = [spmdir '/templates/T1.nii'];
        case 'SPM12'
            tmpcfg.template = [spmdir '/toolbox/OldNorm/T1.nii'];
            addpath([spmdir '/toolbox/OldNorm']);
        otherwise
            error('unknown SPM version')
    end
    
    
    spm_figure('Create','Graphics');
    spm_figure('Redraw');
%     tmpcfg.write = 'yes';
%     tmpcfg.name = normmripath;
%     normalise = ft_volumenormalise(tmpcfg, mri);
    
    templateweight = [];  % optional mask
    mriweight = []; % <---- optional mask, useful for abnormal MRIs
    
    snname = ['s' cfgnemo.participant '_sn.mat'];
    
    params = spm_normalise(tmpcfg.template,...
        mripath, snname,templateweight, mriweight);
    rflags = [];
    rflags.bb = [-85 -120 -60; 85 85 95]; % adjust bounding box, if desired
    
    spm_write_sn(mripath, params, rflags);
end
    
