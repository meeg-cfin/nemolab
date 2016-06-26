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
baselinewindow = [-0.225 -0.025];  % <====== MUST be customized for particular experiment!! ****************
% first band is for ERF
%freqbands = [1 145; 8 13; 13 30; 30 45; 55 75; 75 95; 105 130; 130 145; 155 195];
freqbands = [30 45; 55 75; 75 95; 105 130; 130 145];
tfstats = 0; % NB: tfstats takes a lot of resources!!

cfgnemo.segmethod = 'ftvolseg';
cfgnemo.gridresolution = 10;
cfgnemo.voxelgridtype = 'mni-ft';
cfgnemo.headmodelstrategy =  'openmeeg';

load('standard_sourcemodel3d10mm'); % loads in sourcemodel (i.e., MNI voxel grid)
cfgnemo.sourcemodel = ft_convert_units(sourcemodel,'mm');


cfgnemo.numlayers = 3; % 3 for 3-layer, 4 for 4-layer
cfgnemo.plotvol = 1; % plot surfaces with sensor positions as a check
saveRAM = true; % try to delete mega-matrices after they're no longer needed



%% load and resample data
[~,cwd] = fileparts(pwd); 
switch(cwd)
    case 'nmt_hilbert_test_apr2016'
        load flashes_eyeopen
        cfgnemo.megchans = {'MEG' '-A45' '-A146' '-A147'};
        ergchan = 'E37';
    case 'files' % Aarhus Elekta test
        load data
        
        cfgnemo.megchans = {'MEG'};
        ergchan = 'EOG001';
        
        cfgnemo.megchans = 13:318;
        ergchan = 1;
end
%megergchans = {cfgnemo.megchans{:} ergchan};
megergchans = [ergchan cfgnemo.megchans];

% older saved data lists the chanunit as 'unknown' for reference channels,
% which breaks some FT functions; replace them with 'T'
data.grad.chanunit(find(strcmp(data.grad.chanunit,'unknown')))={'T'}

cfg = [];
cfg.resamplefs = [576];
data = ft_resampledata(cfg, data);

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
%cfgnemo.headmodelstrategy = 'bemcp';
cfgnemo.headmodelstrategy = 'openmeeg';
[leadgrid,vol] = nemo_makeleadfield(cfgnemo);



%% run through filter bank, and obtain timelock and hilbert transform
parfor ii=1:size(freqbands,1)
    cfg=[];
    cfg.channel = megergchans;
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
    cfg.bpfilter = 'yes'; % NB: default butterworth for quick testing; specify more advanced filter for real analysis!
    cfg.bpfreq = freqbands(ii,:);
    cfg.hilbert = 'complex';
    datahilb{ii} = ft_preprocessing(cfg,data);
    
    % filtering first then snipping to shorter time interval avoids edge artifacts
    cfg = [];
    cfg.toilim = toilim;
    datahilb{ii} = ft_redefinetrial(cfg,datahilb{ii});

    databp = datahilb{ii};
    for jj=1:length(databp.trial)
        databp.trial{jj} = real(databp.trial{jj});
    end
    
    

    cfgtl                  = [];
    cfgtl.covariance       = 'yes';
    cfgtl.covariancewindow = 'all';  % may need to change if not desirable to include pre-stim interval
    cfgtl.vartrllength     = 2;
    timelockbp{ii}           = ft_timelockanalysis(cfgtl, databp);

    cfg = [];
    cfg.channel = ergchan;
    erghilb{ii} = ft_preprocessing(cfg,datahilb{ii});
end
if(saveRAM)
    clear data databp
end


%% create spatial filter using the non-Hilbert data
parfor ii=1:size(freqbands,1)
    cfg                   = [];
    cfg.channel           = cfgnemo.megchans;
    cfg.method            = 'lcmv';
    cfg.grid              = leadgrid; % leadfield, which has the grid information
    cfg.vol               = vol; % volume conduction model (headmodel) <-- FIXME: ft_sourceanalysis insists on this even if not necessary (i.e., grid already computed)
    cfg.keepfilter        = 'yes';
    cfg.lcmv.reducerank   = 'no';
    cfg.lcmv.fixedori     = 'yes'; % project on axis of most variance using SVD
    cfg.lcmv.projectnoise = 'yes';
    cfg.lcmv.weightnorm   = 'nai'; %% NOTE: this is crucial for good performance!!! (noise-scaled version of old 'weightnorm')
    cfg.lcmv.keepfilter   = 'yes';
    %    cfg.lcmv.lambda      = '10%';
    source_bp{ii}         = ft_sourceanalysis(cfg, timelockbp{ii}); % high-frequencies only
end

%% apply this spatial filter now to the complex-valued hilbert transform of the data
disp(['Hilbert-sourcing started at ' datestr(now)])
tstart = tic;
for ii=1:size(freqbands,1)
    % change to ordinary "for" if parfor claims "Attempt to serialize data which is too large."
    source_hilbtmp{ii} = ft_apply_spatialfilter(datahilb{ii},source_bp{ii});
    
    datahilb{ii} = []; % don't use "clear" since this is a for loop!!
    
    % perform stats on Hilbert source trials
    nemo_hilbertstats
end
disp(['Hilbert-sourcing finished at ' datestr(now) ' and took ' num2str(toc(tstart)/60) ' minutes to run']);

if(saveRAM)
    clear datahilb
end

%%
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



%% assembles source_hilb{:} into a composite source_tf structure
source_tf=source_hilb{1};
source_tf.freqbands = freqbands(1:end,:);
%source_tf.freqbands(1,:) = [0 0]; % first "band" is actually ERF, this tells the viewer to handle it properly
source_tf.freq = mean(source_tf.freqbands'); % TODO: is this needed for anything??
for jj = 1:size(freqbands,1)
    for ii=1:length(inside_idx)
        source_tf.avg.mom{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.mom{inside_idx(ii)};
        source_tf.avg.itc{inside_idx(ii)}(jj,:)=source_hilb{jj}.avg.itc{inside_idx(ii)};
        if(tfstats)
            source_tf.stat{inside_idx(ii)}(jj,:)=source_hilb{jj}.stat{inside_idx(ii)};
            source_tf.pval{inside_idx(ii)}(jj,:)=source_hilb{jj}.pval{inside_idx(ii)};
        end
    end
end
source_tf.pos = cfgnemo.sourcemodel.pos; % supply MNI pos
source_tf.coordsys = 'mni';

save source_tf source_tf
%%
nemo_plotsourcetf