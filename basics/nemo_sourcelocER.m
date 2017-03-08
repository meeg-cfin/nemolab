% user-defined options:
% data/MRI location
% freqbands: frequency bands of interest
% cfgnemo.segmethod: 'ftvolseg' or 'spm8newseg' or ???
% voxelgridtype: 'mni-ft' is the only smart option :-)
% cfgnemo.headmodelstrategy: 'singleshell', 'openmeeg', 'dipoli', etc.
% cfgnemo.numlayers: number of layers for MRI segmentation (3 or 4 is typical)
%
%
% NOTES
% **** In contrast to 'typical' FieldTrip pipelines, the MEG/EEG is coregistered to the MRI, rather than vice versa
% * Below, a simple Butterworth bandpass is applied, for a quick first look
%   Final results should use sequential lowpass/highpass FIR filters!!
%
% TODO
% - allow choice of reconstruction methods with simple switch (i.e., lcmv vs sloreta/dspm)
% - implement FilterM method
% - add transparent compatibility with EEG (i.e., pass 'elec' vs 'grad' where appropriate)
% - for testing: create test source structures to ensure that nutmegtrip doesn't break for a particular type of data
%                including: lead field plotting, topography plotting, etc.

%%
nemo_ftsetup  % add paths etc.

%% user-defined parameters
toilim = [-0.225 0.25];  % <====== MUST be customized for particular experiment!! ****************
noiselim = [-0.75 -0.3];  % <====== MUST be customized for particular experiment!! ****************
baselinewindow = [-0.225 -0.005];  % <====== MUST be customized for particular experiment!! ****************
activewindow = [0.005 0.225];  % <====== MUST be customized for particular experiment!! ****************
% first band can be wideband; e.g., for ERF
cfgnemo.sourcemethod = 'lcmv';
cfgnemo.headmodelstrategy =  'openmeeg';

cfgnemo.participant = 'SSD';
cfgnemo.stats = 0;
cfgnemo.statori=0;
cfgnemo.VOeyes = 0;
cfgnemo.segmethod = 'ftvolseg';
cfgnemo.numlayers = 3; % 3 for 3-layer, 4 for 4-layer
cfgnemo.plotvol = 0; % plot surfaces with sensor positions as a check
cfgnemo.gridmethod = 'MNI';

load('standard_sourcemodel3d10mm'); % loads in sourcemodel (i.e., MNI voxel grid)
cfgnemo.sourcemodel = ft_convert_units(sourcemodel,'mm');


%%
cfgnemo.megchans = {'MEG'};
megergchans = {cfgnemo.megchans{:}};
%megergchans = {'EMG002', cfgnemo.megchans{:}};

%%
saveRAM = false; % try to delete mega-matrices after they're no longer needed

load('standard_sourcemodel3d10mm'); % loads in sourcemodel (i.e., MNI voxel grid)
sourcemodel = ft_convert_units(sourcemodel,'mm');


%% load and resample data
load flash_1ms_botheyes.mat

% here we filter data that is already segmented; however, it is probably
% desirable to filter the entire dataset first and then segment, to minimize edge artifacts

%%
nemo_mriproc
%%
grad = data.grad;  % NOTE data.hdr.grad does not contain correct information if synthetic gradient has been manipulated above
grad_mm = ft_convert_units(grad,'mm'); % transforms the grad units (m) to the same than mri (mm)
cfgnemo.grad_mri = ft_transform_sens(coreg.meg2mri_tfm, grad_mm); % transforms the grad coordinates to mri coordinates
cfgnemo.grad_mri.coordsys = 'spm';

%%
cfgnemo.bnd = bnd;
if(exist('leadgrid','file'))
    load leadgrid
else
    [leadgrid,vol] = nemo_makeleadfield(cfgnemo);
end


%%
  cfg=[];
    cfg.channel = megergchans;
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
        cfg.hpfilter = 'yes'; % NB: default butterworth for quick testing; specify more advanced filter for real analysis!
        cfg.hpfreq = [75];
        cfg.hpfreq = [1];
        cfg.dftfilter = 'yes';
        cfg.dftfrq = [50:50:200];
    databp = ft_preprocessing(cfg,data);
    
    % filtering first then snipping to shorter time interval avoids edge artifacts
    cfg = [];
    cfg.toilim = toilim;
    databp = ft_redefinetrial(cfg,databp);
    
    
    if(strcmp(cfgnemo.megchans,'MEG')) % correct for magnetometer vs gradiometer scaling
        trials = cell2mat(databp.trial);
        smag = svd(cov(trials(1:3:end,:)'));
        sgrad = svd(cov(trials([2:3:end 3:3:end],:)'));
        
        scaling = sqrt(sgrad(end)/smag(end));
        for ii=1:length(databp.trial)
            databp.trial{ii}(1:3:end,:)=scaling*databp.trial{ii}(1:3:end,:);
        end
        
        inside_idx=find(leadgrid.inside);
        for ii=1:length(inside_idx)
            leadgrid.leadfield{inside_idx(ii)}(1:3:end,:)=scaling*leadgrid.leadfield{inside_idx(ii)}(1:3:end,:);
        end
    end
    
    cfgtl                  = [];
    cfgtl.covariance       = 'yes';
    cfgtl.covariancewindow = 'all';  % may need to change if not desirable to include pre-stim interval
    cfgtl.vartrllength     = 2;
    cfgtl.keeptrials       = 'no';
    timelockbp           = ft_timelockanalysis(cfgtl, databp);
  
%%
    dat=[];
    for jj=1:length(databp.trial)
        dat.b(:,:,jj) = databp.trial{jj};
    end
    Nsamples = size(dat.b,2);
    dat.b = reshape(dat.b,size(dat.b,1),size(dat.b,2)*size(dat.b,3));

    timelockbp.covorig=timelockbp.cov;
%    timelockbp.covfmcd = robustcov(dat.b','Method','fmcd');
%    timelockbp.covogk = robustcov(dat.b','Method','ogk');
%    timelockbp.covoh = robustcov(dat.b','Method','olivehawkins');
%timelockbp.cov = timelockbp.covogk;
    clear dat
    if(0)  % covariance experiments
        timelockbp.covorig = timelockbp.cov;
        Rn = tlbl.cov;
        timelockbp.cov = Rn^(-0.5) * timelockbp.covorig * Rn^(-0.5);
    end

%% create spatial filter using the non-Hilbert data

source_bp=[];
lambda = 10.^(3)
for ii=1:length(lambda)
    cfg                   = [];
    cfg.channel           = {'MEG'};
    cfg.grid              = leadgrid; % leadfield, which has the grid information
    cfg.vol               = vol; % volume conduction model (headmodel) <-- FIXME: ft_sourceanalysis insists on this even if not necessary (i.e., grid already computed)
    cfg.keepfilter        = 'yes';
    cfg.method            = cfgnemo.sourcemethod;
    if strcmp(cfg.method,'sloreta')
        cfg.(cfg.method).lambda = [num2str(lambda(ii)) '%']; %'0%';
    end
     cfg.(cfg.method).reducerank   = 'no';
     if(cfgnemo.statori)
         cfg.grid.mom = zeros(size(cfg.grid.pos))';
         cfg.grid.mom(:,cfg.grid.inside) = ori{ii}';
         cfg.(cfg.method).fixedori     = 'no'; % project on axis of most variance using SVD
     else
         cfg.(cfg.method).fixedori = 'yes';
     end
     cfg.(cfg.method).projectnoise = 'yes';
     cfg.(cfg.method).weightnorm   = 'nai';
    cfg.(cfg.method).keepfilter   = 'yes';
    cfg.keeptrials = 'yes';
    source_bp{ii}         = ft_sourceanalysis(cfg, timelockbp); % high-frequencies only

    % try to make polarities more consistent
    source_bp{ii} = nmt_polaritytweak(source_bp{ii});

    
    inside_idx = find(source_bp{ii}.inside);
    if(cfgnemo.stats)
        clear trials
        source_trials = ft_apply_spatialfilter(databp,source_bp{1});
        for jj=1:length(source_trials.trial)
            trials(jj,:,:) = cell2mat(source_trials.trial(jj).mom);
        end
        clear source_trials

        baselinewindow_idx = dsearchn(source_bp{ii}.time',baselinewindow');

        source_bp{ii}.stat = source_bp{ii}.avg.mom;
        source_bp{ii}.pval = source_bp{ii}.avg.mom;
         
        for jj=1:size(trials,2) % loop over voxels
            voxbl = trials(:,jj,baselinewindow_idx(1):baselinewindow_idx(2));
            voxbl = voxbl(:);
            parfor kk=1:size(trials,3) % loop over time samples
                [p(jj,kk),~,stats]=ranksum(trials(:,jj,kk),voxbl);
                z(jj,kk)=stats.zval;
            end
            source_bp{ii}.stat{inside_idx(jj)}=z(jj,:);
            source_bp{ii}.pval{inside_idx(jj)}=log10(p(jj,:));
        end
        
        
         
    end
    
    
    source_bp{ii}.pos = sourcemodel.pos;
    source_bp{ii}.coordsys = 'spm';
end

src = cell2mat(source_bp{1}.avg.mom);
max(abs(src(:)))

%%
if(1)
%for ii=1:length(lambda)
    cfgplot = [];
    cfgplot.funparameter = 'mom';
    cfgplot.mripath = 'wsSSD.nii';
    cfgplot.atlas = ft_read_atlas('ROI_MNI_V4.nii'); % AAL atlas
    nmt_sourceplot(cfgplot,source_bp{ii});
%    pause
%end
end
