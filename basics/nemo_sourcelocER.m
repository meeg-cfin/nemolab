% NEMO Evoked Response Source Localization pipeline template
%
% NOTES
% **** In contrast to 'typical' FieldTrip pipelines, the MEG/EEG is
%      coregistered to the MRI, rather than vice versa. This is done so
%      that the same participant has a consistent reference MRI, coordinate
%      space, and voxel coordinates, even if sensor positions change across
%      different datasets. It also allows computationally intensive forward
%      models (e.g. BEM or FEM) to be computed once per participant and
%      simply reused with different sensor positions or VOI selections
% * Below, a simple Butterworth bandpass is applied, for a quick first look
%   Final results should use sequential lowpass/highpass FIR filters!!
%
% TODO
% - allow choice of reconstruction methods with simple switch (i.e., lcmv vs sloreta/dspm)
% - add transparent compatibility with EEG (i.e., pass 'elec' vs 'grad' where appropriate)
% - implement SimBio method, requires tetrahedral mesh option
% - for testing: create test source structures to ensure that nutmegtrip doesn't break for a particular type of data
%                including: lead field plotting, topography plotting, etc.

%%
nemo_ftsetup  % add paths etc.

%% user-defined parameters
cfgnemo.participant = 'SSD';

toilim = [-0.225 0.25];  % <====== MUST be customized for particular experiment!! ****************
noiselim = [-0.75 -0.3];  % <====== MUST be customized for particular experiment!! ****************
baselinewindow = [-0.1 -0.005];  % <====== MUST be customized for particular experiment!! ****************
activewindow = [0.005 0.225];  % <====== MUST be customized for particular experiment!! ****************

cfgnemo.sourcemethod = 'lcmv';
cfgnemo.headmodelstrategy =  'openmeeg'; % typically 'openmeeg' or 'singleshell'
%cfgnemo.headmodelstrategy =  'singleshell';

cfgnemo.stats = 0;
cfgnemo.statori = 0; % uses statistical determination of optimum orientation, instead of canonical formula
cfgnemo.segmethod = 'ftvolseg';
cfgnemo.gridmethod = 'MNI';
cfgnemo.numlayers = 3; % 3 for 3-layer, 4 for 4-layer
cfgnemo.plotvol = 0; % plot surfaces with sensor positions as a check
cfgnemo.VOeyes = 0; % experimental; includes eyes in the "brain" VOI

saveRAM = false; % try to delete mega-matrices after they're no longer needed

%% select channels
cfgnemo.megchans = {'MEGMAG'};
megergchans = {cfgnemo.megchans{:}};
%megergchans = {'EMG002', cfgnemo.megchans{:}};

%% load standard head model in desired voxel size (e.g. 4mm, 5mm, or 10mm)
load('standard_sourcemodel3d4mm'); % loads in sourcemodel (i.e., MNI voxel grid)
cfgnemo.sourcemodel = ft_convert_units(sourcemodel,'mm');

%% load and optionally resample data
load flash_1ms_botheyes.mat

%% segment and mesh MRI surfaces
nemo_mriproc
%%
grad = data.grad;  % NOTE data.hdr.grad does not contain correct information if synthetic gradient has been manipulated above
grad_mm = ft_convert_units(grad,'mm'); % transforms the grad units (m) to the same than mri (mm)
cfgnemo.grad_mri = ft_transform_sens(coreg.meg2mri_tfm, grad_mm); % transforms the grad coordinates to mri coordinates
cfgnemo.grad_mri.coordsys = 'spm';

%%
cfgnemo.bnd = bnd;
if(exist([cfgnemo.participant '_leadgrid.mat'],'file'))
    load([cfgnemo.participant '_leadgrid.mat'])
    load([cfgnemo.participant '_headmodel.mat'])
else
    [leadgrid,headmodel] = nemo_makeleadfield(cfgnemo);
end


%% pre-process data (e.g. filtering, segmenting)
  cfg=[];
%    cfg.channel = megergchans;
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
    cfg.hpfilter = 'yes'; % NB: default butterworth for quick testing; specify more advanced filter for real analysis!
    cfg.hpfreq = [1];
    cfg.dftfilter = 'yes'; % powerline filter if needed
    cfg.dftfrq = [50:50:200];
    databp = ft_preprocessing(cfg,data);
    
    % filtering first then snipping to shorter time interval avoids edge artifacts
    cfg = [];
    cfg.toilim = toilim;
    databp = ft_redefinetrial(cfg,databp);
    
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

% try alternative "robust" covariance calculations
%    timelockbp.covorig=timelockbp.cov;
%    timelockbp.covfmcd = robustcov(dat.b','Method','fmcd');
%    timelockbp.covogk = robustcov(dat.b','Method','ogk');
%    timelockbp.covoh = robustcov(dat.b','Method','olivehawkins');
%    timelockbp.cov = timelockbp.covogk;
    clear dat
%     if(0)  % use actual baseline to estimate noise covariance
%         timelockbp.covorig = timelockbp.cov;
%         Rn = tlbl.cov;
%         timelockbp.cov = Rn^(-0.5) * timelockbp.covorig * Rn^(-0.5);
%     end

%% create spatial filter using the non-Hilbert data

if(0)
    warndlg('reducerank = 2 !!!')
    cfgred = [];
    cfgred.reducerank = 2;
    leadgridred = nmt_reducerank(cfgred,leadgrid);
else
    cfgred.reducerank = 3;
    leadgridred = leadgrid;
end

source_bp=[];
for ii=1:length(lambda)
    cfg                   = [];
    cfg.channel           = {'MEGMAG'};
    cfg.grid              = leadgridred; % leadfield, which has the grid information
    cfg.headmodel         = headmodel; % volume conduction model (headmodel) <-- FIXME: ft_sourceanalysis insists on this even if not necessary (i.e., grid already computed)
    cfg.keepfilter        = 'yes';
    cfg.method            = cfgnemo.sourcemethod;
    if strcmp(cfg.method,'sloreta')
        lambda = 10.^(8)
        cfg.(cfg.method).lambda = [num2str(lambda(ii)) '%']; %'0%';
    else
        cfg.(cfg.method).projectnoise = 'yes';
        cfg.(cfg.method).weightnorm   = 'nai';
    end
    if(cfgnemo.statori)
        cfg.grid.mom = zeros(size(cfg.grid.pos))';
        cfg.grid.mom(:,cfg.grid.inside) = ori{ii}';
        cfg.(cfg.method).fixedori     = 'yes'; % project on axis of most variance using SVD
    else
        cfg.(cfg.method).fixedori = 'yes'
    end
    cfg.(cfg.method).keepfilter   = 'yes';
    cfg.keeptrials = 'yes';
    source_bp{ii}         = ft_sourceanalysis(cfg, timelockbp); % high-frequencies only

    inside_idx = find(source_bp{ii}.inside);

    % try to make polarities more consistent
    if(strcmp(cfg.(cfg.method).fixedori,'yes'))
        cfgpol = [];
        cfgpol.toilim = [0 0.2];
        source_bp{ii} = nmt_polaritytweak(cfgpol,source_bp{ii});
    else
        source_bp{ii}.avg.mompow = source_bp{ii}.avg.mom;
        source_bp{ii}.avg.ori = source_bp{ii}.avg.mom;
        for jj=1:length(inside_idx)
            source_bp{ii}.avg.mompow{inside_idx(jj)} = sqrt(sum((source_bp{ii}.avg.mom{inside_idx(jj)}).^2));
            source_bp{ii}.avg.ori{inside_idx(jj)} = source_bp{ii}.avg.mom{inside_idx(jj)};
        end
    end
    
    
   
    if(0)
        source_bp{ii}.avg.momhilb = source_bp{ii}.avg.mom;
        source_bp{ii}.avg.test = source_bp{ii}.avg.mom;
        for jj=1:length(inside_idx)
            source_bp{ii}.avg.momhilb{inside_idx(jj)} = abs(hilbert(source_bp{ii}.avg.mom{inside_idx(jj)}));
            noisethresh = std(source_bp{ii}.avg.mom{inside_idx(jj)}(1:220));
            source_bp{ii}.avg.noisethresh{inside_idx(jj)} = noisethresh;
        end
    end
    
    
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
    
    
    source_bp{ii}.pos = cfgnemo.sourcemodel.pos;
    source_bp{ii}.coordsys = 'spm';
end

src = cell2mat(source_bp{1}.avg.mom);
max(abs(src(:)))  % display max amplitude of source output

%% plot results
if(1)
    cfgplot = [];
    
    if(strcmp(cfg.(cfg.method).fixedori,'yes'))
        cfgplot.funparameter = 'mom';
    else
        cfgplot.funparameter = 'mompow';
    end
    %    cfgplot.vecparameter = 'ori';
    cfgplot.mripath = 'wsSSD.nii';
    cfgplot.atlas = ft_read_atlas('ROI_MNI_V4.nii'); % AAL atlas
    cfgplot.colormap = cmocean('balance',64); % requires cmocean colormap toolbox
    
    if(1) % masking, e.g. for thresholding source maps
        cfgplot.maskparameter = 'msk';
        cfgplot.masktype = 'momabs';  % masktypes: pval, momhigh, momlow, momabs, itchigh, itclow, pow
        maskthresh = 3; % for pval masktype, this should be log(pval)
        
        switch(cfgplot.funparameter) % works for pow but not mom...
            case {'mom','momhilb'} % FIXME: this is really messy, because of the cell strucutre of avg.mom :-(
                % nmt_polaritytweak;
                mom = reshape([source_bp{ii}.avg.mom{:}],length(source_bp{ii}.time),length(inside_idx));
                [peakval,peakvox]=max(abs(mom),[],2);
                
            case 'mompow'
                source_bp{ii}.msk = zeros(size(source_bp{ii}.avg.mom,1),1,length(source_bp{ii}.time));
                for tt=1:length(source_bp{ii}.time) % pick threshold relative to peak at each given time point
                    source_bp{ii}.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*max(abs(mom(:)));
                end
        end
        
        source_bp{ii}.msk = zeros(size(source_bp{ii}.avg.mom,1),1,length(source_bp{ii}.time));
        switch(cfgplot.masktype)
            case 'noisethresh'
                noisethresh=repmat([source_bp{ii}.avg.noisethresh{:}],size(mom,1),1);
                for jj=1:size(mom,1)
                    momsort = sort(abs(mom(jj,:)),'descend');
                    momthresh(jj,:) = abs(mom(jj,:)) > momsort(round(.02*(size(mom,2))));
                end
                source_bp{ii}.msk(inside_idx,1,:) = (abs(mom') > 3.5*noisethresh') & momthresh';
            case 'momhigh'
                source_bp{ii}.msk(inside_idx,1,:) = mom' > maskthresh;
            case 'momlow'
                source_bp{ii}.msk(inside_idx,1,:) = mom' < maskthresh;
            case 'momabs'
                source_bp{ii}.msk(inside_idx,1,:) = abs(mom') > maskthresh;
            case 'momratio'
                for tt=1:length(source_bp{ii}.time) % pick threshold relative to peak at each given time point
                    %            source_bp{ii}.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*peakval(tt);
                    source_bp{ii}.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*max(abs(mom(:)));
                end
            otherwise
                error('masktype unknown');
        end
    end
    
    nmt_sourceplot(cfgplot,source_bp{ii});
end
