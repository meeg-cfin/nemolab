% user-defined options:
% data/MRI location
% freqbands: frequency bands of interest
% cfgnemo.segmethod: 'ftvolseg' or 'spm8newseg' or ???
% voxelgridtype: 'mni-ft' is the only smart option :-)
% tfstats: hilbert stats on/off (can require lots of RAM!!)
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
% - investigate whether cfg.precision = 'single' can save memory without causing problems!
% - *** CHECK ITC COMPUTATION
% - *** HOW TO HANDLE ERF COMPUTATION?? (does running it through hilbert make sense??)
% - allow choice of reconstruction methods with simple switch (i.e., lcmv vs sloreta/dspm)
% - allow lead field choice
% - implement SimBio method, requires tetrahedral mesh option
% - downsampling to save on memory/computation
% - check that large matrices are cleared when they are no longer needed
% - implement FilterM method
% - compare with dipoli
% - add transparent compatibility with EEG (i.e., pass 'elec' vs 'grad' where appropriate)
% - allow manual SPM normalization, and save result
% - add "block design" mode and nutmegtrip support for it
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
%freqbands = [1 145; 8 13; 13 30; 30 45; 55 75; 75 95; 105 130; 130 145; 155 195];
%freqbands = [1 145; 30 45; 55 75; 75 95; 105 130; 130 145; 155 195; 205 245; 255 295];
freqbands = [55 75; 75 95; 105 130; 130 145];
%freqbands = [55 75; 75 95; 105 120; 120 145; 155 195; 205 245; 255 295; 305 345; 355 395];
%freqbands = [1 130];
cfgnemo.tfstats = 1; % NB: tfstats may take lots of resources!!
cfgnemo.statori = 0;
cfgnemo.covmethod = 1; % 0 = 'ordinary' cov; 1 = median cov; 2 = 'robust' cov
cfgnemo.downsample = 1;
cfgnemo.segmethod = 'ftvolseg';
cfgnemo.voxelgridtype = 'mni-ft';
cfgnemo.sourcemethod = 'lcmv';
cfgnemo.headmodelstrategy =  'openmeeg'; %'simbio';
cfgnemo.VOeyes = 0; % include eyes in VOI
cfgnemo.sourceplottype = 'fancy';

load('standard_sourcemodel3d10mm'); % loads in sourcemodel (i.e., MNI voxel grid)
cfgnemo.sourcemodel = ft_convert_units(sourcemodel,'mm');


cfgnemo.numlayers = 3; % 3 for 3-layer, 4 for 4-layer
cfgnemo.plotvol = 0; % plot surfaces with sensor positions as a check
saveRAM = false; % try to delete mega-matrices after they're no longer needed


%% load and resample data
[~,cwd,ext] = fileparts(pwd); 
cwd = [cwd ext];
switch(cwd)
    case 'nmt_hilbert_test_apr2016'
        load flashes_eyeopen
        badchans = {'-A45' '-A146' '-A147'};
        cfgnemo.megchans = {'MEG' badchans{:}};
        cfgnemo.ergchan = 'E37';
        cfgnemo.participant = 'SSD';
    case 'flashesL_CM43403'
        load data
        cfgnemo.megchans = {'MEG'};
        cfgnemo.ergchan = 'E37';
        cfgnemo.participant = 'CM';
    case '001.flash_01_raw'
        load data
        badchans = {'-MEG2511'};
        cfgnemo.ergchan = 'EOG001';
        
        cfgnemo.megchans = {'MEGMAG' badchans{:}};
        cfgnemo.participant = 'SSD';
    case {'001.flash_1ms_righteye','002.flash_1ms_lefteye','003.flash_1ms_botheyes','001.flash_1ms_botheyes','002.flash_5ms_botheyes'}
        cfg=[];
        datafile = dir('*.fif');
        cfg.dataset = datafile.name;
        cfg.channel = 'MISC001';
        cfg.coilaccuracy = 1;
        data=ft_preprocessing(cfg);
        
        [pks,pkidx]=findpeaks(data.trial{1},'MinPeakProminence',.1); %.022 for MP
        
        badchans = {};
        cfgnemo.megchans = {'MEGMAG' badchans{:}};
        cfg.trl = [pkidx'-2500 pkidx'+2500];
        cfg.trl(:,3) = -2500;
        % in recent datasets, "EMG" is actually "ERG"
        cfg.channel = {'EOG001' 'EOG002' 'EMG001' 'EMG002' cfgnemo.megchans{:}};
        cfg.dftfilter = 'no';
        data = ft_preprocessing(cfg);
        switch(cwd)
            case {'001.flash_1ms_botheyes','002.flash_5ms_botheyes'}
                cfgnemo.ergchan = 'EMG002';
            case '001.flash_1ms_righteye'
                cfgnemo.ergchan = 'EOG001';
            otherwise
                cfgnemo.ergchan = 'EOG002';
        end
        
        
        cfgnemo.participant = 'SSD';

end
megergchans = {cfgnemo.megchans{:} cfgnemo.ergchan};


% older saved data lists the chanunit as 'unknown' for reference channels,
% which breaks some FT functions; replace them with 'T'
data.grad.chanunit(find(strcmp(data.grad.chanunit,'unknown')))={'T'}

if(cfgnemo.downsample) % optionally downsample
    cfg = [];
    cfg.resamplefs = [1000];
    data = ft_resampledata(cfg, data);
end



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
[leadgrid,vol] = nemo_makeleadfield(cfgnemo);

%%
cfg=[];
cfg.channel = megergchans;
cfg.bpfilter = 'yes';
cfg.bpfreq = [3 195];
evoked = ft_timelockanalysis(cfg,data);

%% run through filter bank, and obtain timelock and hilbert transform
parfor ii=1:size(freqbands,1)
    cfg=[];
    cfg.channel = megergchans;
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
    if(1)
        cfg.bpfilter = 'yes'; % NB: default butterworth for quick testing; specify more advanced filter for real analysis!
        cfg.bpfreq = freqbands(ii,:);
    else
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [freqbands(ii,1) 195];
    end
    cfg.hilbert = 'complex';
    datahilb{ii} = ft_preprocessing(cfg,data);
    
    % filtering first then snipping to shorter time interval avoids edge artifacts
    cfg = [];
    cfg.toilim = toilim;
    datahilb{ii} = ft_redefinetrial(cfg,datahilb{ii});
    
    cfg = [];
    cfg.toilim = baselinewindow;
    controlbp = ft_redefinetrial(cfg,datahilb{ii});
    
    cfg = [];
    cfg.toilim = activewindow;
    activebp = ft_redefinetrial(cfg,datahilb{ii});
    
    
    databp = datahilb{ii};
    for jj=1:length(databp.trial)
        databp.trial{jj} = real(databp.trial{jj});
    end
    
    for jj=1:length(controlbp.trial)
        controlbp.trial{jj} = real(controlbp.trial{jj});
    end
    
    for jj=1:length(controlbp.trial)
        activebp.trial{jj} = real(activebp.trial{jj});
    end
    
    
    
    cfgtl                  = [];
    cfgtl.covariance       = 'yes';
    cfgtl.covariancewindow = 'all';  % may need to change if not desirable to include pre-stim interval
    cfgtl.vartrllength     = 2;
%    cfgtl.keeptrials       = 'yes';
    timelockbp{ii}           = ft_timelockanalysis(cfgtl, databp);
    timelockcontrol{ii}           = ft_timelockanalysis(cfgtl, controlbp);
    timelockactive{ii}           = ft_timelockanalysis(cfgtl, activebp);
    
    cfg = [];
    cfg.channel = cfgnemo.ergchan;
    erghilb{ii} = ft_preprocessing(cfg,datahilb{ii});
end
if(saveRAM)
    clear data databp
end

baselinewindow_idx = dsearchn(datahilb{1}.time{1}',baselinewindow');
activewindow_idx = dsearchn(datahilb{1}.time{1}',activewindow');




%% median covariance -- you definitely don't have enough RAM to parfor this!!!!!!!!
if(cfgnemo.covmethod==1)
   
    meglabels=ft_channelselection(cfgnemo.megchans,leadgrid.label);
    [~,ommagid] = intersect(leadgrid.label,meglabels);
    [~,ftmegchanid] = intersect(datahilb{1}.label,meglabels);

    % chop out bits of trial that aren't used for contrast?
    %tsel = [baselinewindow_idx(1):baselinewindow_idx(2) activewindow_idx(1):activewindow_idx(2)];
    
    for ii=1:size(freqbands,1)
        disp(freqbands(ii,:))
        
        dat=[];
        for jj=1:length(datahilb{ii}.trial)
            dat.b(:,:,jj) = real(datahilb{ii}.trial{jj}(ftmegchanid,:));
        end
        Nsamples = size(dat.b,2);
        dat.b = reshape(dat.b,size(dat.b,1),size(dat.b,2)*size(dat.b,3));
%        timelockbp{ii}.cov=zeros(1,size(dat.b,1),size(dat.b,2)); % reinitialize cov matrix (removes original unaveraged covariance)
        timelockbp{ii}.cov(ftmegchanid,ftmegchanid) = robustcov(dat.b','Method','olivehawkins');
        
        % generate covariances for active and control periods, for possible
        % use in orientation selection
        timelockbp{ii}.covact=zeros(size(timelockbp{ii}.cov));
        timelockbp{ii}.covcon=zeros(size(timelockbp{ii}.cov));
        trialstart_idx = [0:Nsamples:length(dat.b)];
        trialstart_idx(end) = []; % last one is too far!
        activesamples = [];
        controlsamples = [];
        for nn=1:length(trialstart_idx)
            activesamples =  [activesamples trialstart_idx(nn)+(activewindow_idx(1):activewindow_idx(2))];
            controlsamples =  [controlsamples trialstart_idx(nn)+(baselinewindow_idx(1):baselinewindow_idx(2))];
        end
        timelockbp{ii}.covcon(ftmegchanid,ftmegchanid) = robustcov(dat.b(:,controlsamples)');
        timelockbp{ii}.covact(ftmegchanid,ftmegchanid) = robustcov(dat.b(:,activesamples)');
        

%                 
%         for jj=1:size(dat.b,2)
%             dat.C(:,:,jj) = (dat.b(:,jj))*(dat.b(:,jj))';
%         end
%         timelockbp{ii}.cov=zeros(1,size(dat.C,1),size(dat.C,2));
%         for jj=1:size(dat.C,1) % to save memory, loop through columns for median
%             timelockbp{ii}.cov(1,:,jj) = median(dat.C(:,jj,:),3);
%         end
%         
%         % generate covariances for active and control periods, for possible
%         % use in orientation selection
%         timelockbp{ii}.covact=zeros(size(timelockbp{ii}.cov));
%         timelockbp{ii}.covcon=zeros(size(timelockbp{ii}.cov));
%         trialstart_idx = [0:Nsamples:length(dat.b)];
%         trialstart_idx(end) = []; % last one is too far!
%         activesamples = [];
%         controlsamples = [];
%         for nn=1:length(trialstart_idx)
%             activesamples =  [activesamples trialstart_idx(nn)+(activewindow_idx(1):activewindow_idx(2))];
%             controlsamples =  [controlsamples trialstart_idx(nn)+(baselinewindow_idx(1):baselinewindow_idx(2))];
%         end
%         for jj=1:size(dat.C,1) % to save memory, loop through columns for median
%             timelockbp{ii}.covcon(1,:,jj) = median(dat.C(:,jj,controlsamples),3);
%         end
%         
%         
        
        %     timelockbp{ii}.cov = reshape(median(dat.C,3),1,size(dat.C,1),size(dat.C,1));
        clear dat
    end
end


%% determine orientation based on direction that maximizes statistic
if(cfgnemo.statori)
    if(exist('./statori.mat','file'))
        load statori
    else
        disp(['SnPM orientation started at ' datestr(now)])
        tstart = tic;
        
        cfg = [];
        cfg.method = 'lcmv';
        cfg.stat = 'powrat';

        parfor ii=1:size(freqbands,1)
            if(~strcmp(cfg.stat,'powrat'))
                error('Aborting to prevent server explosion :-) This is insanity with parfor, change to normal for!');
            end
            disp(freqbands(ii,:))
            inside_idx=find(leadgrid.inside);
            
            meglabels=ft_channelselection(cfgnemo.megchans,leadgrid.label);
            [~,ommagid] = intersect(leadgrid.label,meglabels);
            [~,ftmegchanid] = intersect(datahilb{ii}.label,meglabels);
            
            if(strcmp(cfg.stat,'ranksumactive'))
                dat=[];
                for jj=1:length(datahilb{ii}.trial)
                    dat.b(:,:,jj) = datahilb{ii}.trial{jj}(ftmegchanid,:);
                end
                dat.controlwin = baselinewindow_idx(1):baselinewindow_idx(2);
                dat.activewin = activewindow_idx(1):activewindow_idx(2);;
                
                dat.bactive = reshape(dat.b(:,dat.activewin,:),size(dat.b,1),length(dat.activewin)*size(dat.b,3));
                dat.bcontrol = reshape(dat.b(:,dat.controlwin,:),size(dat.b,1),length(dat.controlwin)*size(dat.b,3));
            end
            
            
            dat.C = squeeze(mean(timelockbp{ii}.cov(:,ftmegchanid,ftmegchanid),1));
            dat.invC = inv(dat.C);

            dat.Cact = squeeze(mean(timelockbp{ii}.covact(:,ftmegchanid,ftmegchanid),1));
            dat.Ccon = squeeze(mean(timelockbp{ii}.covcon(:,ftmegchanid,ftmegchanid),1));
            dat.invCact = inv(dat.Cact);
            dat.invCcon = inv(dat.Ccon);
            
            
            
            ft_progress('init','etf');
            for kk=1:length(inside_idx)
                dat.L = leadgrid.leadfield{inside_idx(kk)}(ommagid,:);
                
                [ori{ii}(kk,:)] = nemo_spatfilt_statori(cfg,dat);
                ft_progress(kk/length(inside_idx),'%d of %d',kk,length(inside_idx));
            end
            ft_progress('close');
        end
        
        save statori ori
        
        disp(['SnPM orientation finished at ' datestr(now) ' and took ' num2str(toc(tstart)/60) ' minutes to run']);
    end
    
end

%% create spatial filter using the non-Hilbert data
for ii=1:size(freqbands,1)

    cfg                   = [];
    cfg.channel           = cfgnemo.megchans;
    cfg.grid              = leadgrid; % leadfield, which has the grid information
    cfg.vol               = vol; % volume conduction model (headmodel) <-- FIXME: ft_sourceanalysis insists on this even if not necessary (i.e., grid already computed)
    cfg.keepfilter        = 'yes';
    cfg.method            = cfgnemo.sourcemethod;
    if strcmp(cfg.method,'sloreta')
        cfg.(cfg.method).lambda = '1000%';
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
     cfg.(cfg.method).weightnorm   = 'nai'; %% NOTE: nai or lfnorm seems crucial for good performance!
     %       cfg.(cfg.method).invCact = inv(squeeze(mean(timelockactive{ii}.cov(:,2:end,2:end),1)));
 %        cfg.(cfg.method).invCcon = inv(squeeze(mean(timelockcontrol{ii}.cov(:,2:end,2:end),1)));;
%    cfg.(cfg.method).weightnorm   = 'naiori'; %% NOTE: nai or lfnorm seems crucial for good performance!
%    cfg.(cfg.method).weightnorm   = 'lfnorm'; %% NOTE: nai or lfnorm seems crucial for good performance!
%cfg.(cfg.method).lambda      = '10%';
    cfg.(cfg.method).keepfilter   = 'yes';
    source_bp{ii}         = ft_sourceanalysis(cfg, timelockbp{ii}); % high-frequencies only
end

%% apply this spatial filter now to the complex-valued hilbert transform of the data
disp(['Hilbert-sourcing started at ' datestr(now)])
tstart = tic;
for ii=1:size(freqbands,1)
    % change to ordinary "for" if parfor claims "Attempt to serialize data which is too large."
    source_hilbtmp{ii} = ft_apply_spatialfilter(datahilb{ii},source_bp{ii});
    
    if(saveRAM)
        datahilb{ii} = []; % don't use "clear" since this is a for loop!!
    end
    
    % perform stats on Hilbert source trials
    nemo_hilbertstats
end
disp(['Hilbert-sourcing finished at ' datestr(now) ' and took ' num2str(toc(tstart)/60) ' minutes to run']);

if(saveRAM)
    clear datahilb
end


%% assembles source_hilb{:} into a composite source_tf structure
source_tf=source_hilb{1};
source_tf.freqbands = freqbands(1:end,:);

% to ensure inclusion of "voxel" containing ERG channel
source_tf.inside = source_hilbtmp{1}.inside;

%source_tf.freqbands(1,:) = [0 10]; % first "band" is actually ERF, this tells the viewer to handle it properly
source_tf.freq = mean(source_tf.freqbands'); % TODO: is this needed for anything??
for jj = 1:size(freqbands,1)
    source_hilb{jj}.avg.ori{inside_idx(1)} = [0;0;0]; % dummy ori for ERG channel
%    avgfields = fieldnames(source_hilb{jj}.avg);
    for ii=1:length(inside_idx)
%         for kk=1:length(avgfields) % TODO: don't know how to make this work yet
%            source_tf.avg.(avgfields{kk}){inside_idx(ii)}{jj}=source_hilb{jj}.avg.(avgfields{kk}){inside_idx(ii)};
%         end

        source_tf.avg.mom{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.mom{inside_idx(ii)};
        source_tf.avg.aa{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.aa{inside_idx(ii)};
        if(cfgnemo.statori==0)
            source_tf.avg.ori{inside_idx(ii)}(:,jj)=source_hilb{jj}.avg.ori{inside_idx(ii)};
        end
        %        source_tf.avg.mom{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.mom{inside_idx(ii)}/max(abs(source_hilb{jj}.avg.mom{inside_idx(ii)}))*max(abs(source_hilb{jj}.stat{inside_idx(ii)}));
        %        source_tf.avg.aa{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.aa{inside_idx(ii)}/max(abs(source_hilb{jj}.avg.aa{inside_idx(ii)}))*max(abs(source_hilb{jj}.stat{inside_idx(ii)}));
        source_tf.avg.itc{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.itc{inside_idx(ii)};
        if(isfield(source_tf.avg,'covcond'))
            source_tf.avg.covcond(jj)=source_hilb{jj}.avg.covcond;
        end
            
        if(cfgnemo.tfstats)
            source_tf.stat{inside_idx(ii)}(jj,:)=source_hilb{jj}.stat{inside_idx(ii)};
            source_tf.pval{inside_idx(ii)}(jj,:)=source_hilb{jj}.pval{inside_idx(ii)};
            source_tf.statitc{inside_idx(ii)}(jj,:)=source_hilb{jj}.statitc{inside_idx(ii)};
            source_tf.pitc{inside_idx(ii)}(jj,:)=source_hilb{jj}.pitc{inside_idx(ii)};
        end
    end
end
source_tf.pos = cfgnemo.sourcemodel.pos; % supply MNI pos
source_tf.coordsys = 'mni';

save source_tf source_tf
%%
nemo_plotsourcetf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% alternative code after this point %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0) % alternative method: this is probably slower
cfgtl                  = [];
cfgtl.channel          = channel;
cfgtl.vartrllength     = 2;
cfgtl.keeptrials       = 'yes';
timelockhilb{ii}       = ft_timelockanalysis(cfgtl, datahilb);
timelockhilb{ii}.covariance = timelockbp{ii}.covariance;

cfg                   = [];
cfg.channel           = channel;
cfg.method            = 'lcmv';
cfg.grid              = leadgrid; % leadfield, which has the grid information
cfg.vol               = vol; % volume conduction model (headmodel)
cfg.keepfilter        = 'yes';
cfg.lcmv.reducerank   = 'no';
cfg.lcmv.fixedori     = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.weightnorm   = 'nai'; %% NOTE: this is crucial for good performance!!! (noise-scaled version of old 'weightnorm')
cfg.lcmv.keepfilter   = 'yes';
%    cfg.lcmv.lambda     = '10%';
source_hilbtmp{ii}    = ft_sourceanalysis(cfg, timelockhilb{ii}); % high-frequencies only
end

%%
if(0) % real part of this same as source_bp{ii}; possibly useful for average phase or some other reason???
    cfg.rawtrial='no';
    source_hilbtmp2 = ft_sourceanalysis(cfg, timelockhilb);
end
