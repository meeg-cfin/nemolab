%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Hilbert and statistics %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variable space needs to have loaded:
% data (renamed below!)

% output of this script can either be plotted with nemo_plot_hilbertchannel (make sure to specify the frequency bands as lower and upper limits per band!)
% or with fieldtrip's own plotting functions, ft_singleplotTFR and ft_multiplotTFR (if topo-data) - then make sure to specify data.freq containing the mean of each 
% frequency band. 

freqbands = [55 75; 75 95; 105 125; 125 145];

% specs for firws filter:
tswidth = 6;

baselinewindow = [-.15 -.02];
toilim  = [-.15 .25]

stim = 'off';
filtertype = 'firws'; % 'but' or 'firws'


% specs for firws filter:
tswidth = 6;

if(strcmp(stim, 'on'))
    data = dataon_clean;              % data already loaded - adjust for your needs. This is preprocessed + downsampled sensorlevel data.
    display('Analyzing ON data');
elseif(strcmp(stim, 'off'))
    data = dataoff_clean;
    display('Analyzing OFF data');
else
    error('unknown stimulus condition!')
end

saving = 1;     % should data be saved?
plotting = 0;   % should results be plotted?

%% run through filter bank, and obtain timelock and hilbert transform


% preallocation
timelockbp = cell(1,size(freqbands,1));
data_hilb_bp  = cell(1,size(freqbands,1));
filterdefinitions = cell(1,size(freqbands,1));

for ii=1:size(freqbands,1)   % using parfor here has some weird effects on filtering feedback etc
    cfg=[];
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
    if(strcmp(filtertype, 'but'))
        cfg.bpfilter = 'yes'; % NB: default butterworth for quick testing; specify more advanced filter for real analysis!
        cfg.bpfreq = freqbands(ii,:);
    elseif(strcmp(filtertype, 'firws'))
        cfg.hpfilter = 'yes';
        cfg.hpfreq = freqbands(ii,1);
        cfg.hpfilttype = filtertype;
        cfg.hpfiltdf = tswidth;
        cfg.hpfiltdir = 'onepass-zerophase';
        cfg.plotfiltresp = 'no';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = freqbands(ii,2);
        cfg.lpfilttype = filtertype;
        cfg.lpfiltdf = tswidth;
        cfg.lpfiltdir = 'onepass-zerophase';
        %         cfg.bpfilter = 'yes';
        %         cfg.bpfreq = freqbands(ii,:);
        %         cfg.bpfilttype = filtertype;
        %         cfg.bpfiltdf = tswidth;
        %         cfg.bpfiltdir = 'onepass-zerophase';
        %         cfg.plotfiltresp = 'yes';
    else
        error('I do not know the specified filter type, please add!')
    end
    databp = ft_preprocessing(cfg,data);
    
    %     if(isfield(databp.cfg, 'printonce'))
    %         filterdefinitions{ii}.low=databp.cfg.printonce.identifier.ft_preprocessing.preproc.ft_preproc_lowpassfilter;   % keep that because databp will be deleted.
    %         filterdefinitions{ii}.high=databp.cfg.printonce.identifier.ft_preprocessing.preproc.ft_preproc_highpassfilter;   % keep that because databp will be deleted.
    %     end
    
    % filtering first then snipping to shorter time interval avoids edge artifacts
    cfg = [];
    cfg.toilim = toilim;
    databp = ft_redefinetrial(cfg,databp);
    
    cfgtl                  = [];
    cfgtl.covariance       = 'yes';
    cfgtl.covariancewindow = 'all';  % may need to change if not desirable to include pre-stim interval
    cfgtl.vartrllength     = 2;
    timelockbp{ii}           = ft_timelockanalysis(cfgtl, databp);
    
    cfghilb=[];
    cfghilb.hilbert = 'complex';
    data_hilb_bp{ii} = ft_preprocessing(cfghilb,databp);
    
    %     cfg = [];
    %     cfg.channel = ergchan;
    %     erghilb{ii} = ft_preprocessing(cfg,datahilb{ii});
end
clear databp


%% apply this spatial filter now to the complex-valued hilbert transform of the data
addpath /Data2/AarhusOnOff/scripts_britta/
addpath /home/britta/Documents/Matlab_toolboxes/circstat/

saveRAM = 0;
tfstats = 1;

disp(['Hilbert-sourcing started at ' datestr(now)])
tstart = tic;
for ii=1:size(freqbands,1)
    % perform stats on sensorlevel
    nemo_hilbertstats_sensorlevel
end
disp(['Hilbert-sourcing finished at ' datestr(now) ' and took ' num2str(toc(tstart)/60) ' minutes to run']);

if(saveRAM)
    clear data_hilb_bp
end


%% assembles source_hilb{:} into a composite source_tf structure
data_tf=data_hilbert{1};
data_tf.avg = zeros(numel(data_hilbert{1}.label), size(freqbands,1), size(data_hilbert{1}.time, 2));
data_tf.itc = data_tf.avg;
data_tf.stat = data_tf.avg;
data_tf.pval = data_tf.avg;


data_tf.freqbands = freqbands;

data_tf.freq = mean(data_tf.freq'); % This is needed if output should be plotted with fieldtrip's plotting functions - ft_singleplotTFR or ft_multiplotTFR
for jj = 1:size(freqbands,1)
    data_tf.avg(:,jj,:)=reshape(data_hilbert{jj}.avg, numel(data_hilbert{1}.label), 1, size(data_hilbert{1}.time, 2));
    data_tf.itc(:,jj,:)=reshape(data_hilbert{jj}.itc, numel(data_hilbert{1}.label), 1, size(data_hilbert{1}.time, 2));
    if(tfstats)
        data_tf.stat(:,jj,:)=reshape(data_hilbert{jj}.stat, numel(data_hilbert{1}.label), 1, size(data_hilbert{1}.time, 2));
        data_tf.pval(:,jj,:)=reshape(data_hilbert{jj}.pval, numel(data_hilbert{1}.label), 1, size(data_hilbert{1}.time, 2));
    end
end

data_tf.dimord = 'chan_freq_time';

if(saving)
    
    if(strcmp(stim, 'on'))
        save(fullfile(outdir, 'erg_tf_ON.mat'), 'data_tf');
        display('Saving ON data');
    elseif(strcmp(stim, 'off'))
        save(fullfile(outdir, 'erg_tf_OFF.mat'), 'data_tf');
        display('Saving OFF data');
    end
end


%% plotting with nemo_plot_hilbertchannel

if(plotting)   
    figure;
      
    cfg = [];
    cfg.funparameter = 'stat';

    a=subplot(2,1,1)
    cfg.channel = 'EOG001';
    plot_hilbertchannel(cfg, data_tf);
    set(a, 'fontsize', 18)
    title('ERG right')
    axis tight
    
    b=subplot(2,1,2)
    cfg.channel = 'EOG002';
    plot_hilbertchannel(cfg, data_tf);
    set(b, 'fontsize', 18)
    title('ERG left')
    axis tight

end
