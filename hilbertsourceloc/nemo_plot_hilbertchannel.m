function nemo_plot_hilbertchannel(cfg, data_tf)

% NEMO_PLOT_HILBERTCHANNEL plots the output of nemo_hilberstats_sensorlevel as TFR of a specified parameter and one channel
%
% Use as
%       plot_hilbertchannel(cfg, data_tf)
%
% You should specify the following parameters:
%   cfg.funparameter       = parameter you want to plot, can be 'stat', 'avg', 'itc', 'pval'
%   cfg.channel            = string specifying the channel that should be plotted, as specified in data_tf.label
%   cfg.freqbands          = n_freq x 2 matrix specifying the lower and upper limits of the frequency bands (see below!)
%   cfg.zlim               = limits for color dimension, can be 'maxabs' (default), 'mincenter' (centered at zero, minimum value is limit), 'maxcenter' or [min max]
%   cfg.xlim               = time limits for plotting
%   cfg.masking            = 0 or 1. If 1, masking is carried out, i.e. the binary mask in cfg.maskparameter will be used to only show the areas where the mask is 1
%   cfg.marking            = 0 or 1. If 1, the areas where the binary mask specified in cfg.maskparameter is 1 will be marked by a black rectangle, if the contain more data points than specified in cfg.markingclust
%   cfg.maskparameter      = if masking or marking is desired, this should point to a field in data_tf that contains a binary mask
%   cfg.markingclust       = number that specifies the minimum number of datapoints in a "cluster" (specified by the binary mask), that should be outlined by a rectangle
%   cfg.zeroline           = 0 or 1. If 1, a vertical line is drawn to indicate 0 on the time axis.
%
% This function needs frequency bands instead of center/mean frequencies.
% You can either have that specified directly in your data as
%   data_tf.freqbands
% or in the cfg as cfg.freqbands

% Author: Britta Westner


% check config structure:
if(~isfield(cfg, 'funparameter'))
    error('Please specify a funparameter in your configuration structure!')
elseif(~isfield(cfg, 'channel'))
    error('Please specify a channel in your configuration structure!')
end

% defaults for masking and marking:
if(~isfield(cfg, 'masking'))
    cfg.masking = 0;
end

if(~isfield(cfg, 'marking'))
    cfg.marking = 0;
end

if(~isfield(cfg, 'zeroline'))
    cfg.zeroline = 0;
end


%% prepare data

% Dimord stuff
% check which dimension is channel
datadimord = tokenize(data_tf.dimord, '_');

if(length(datadimord)~=3)
    error('can only handle 3-dimensional data')
end

if(sum(ismember(datadimord, {'chan', 'time', 'freq'}))~=3)
    error(sprintf('Can only handly time-frequency data on channel-level at the moment. Your dimord is %s.', data_tf.dimord));
end


% Channel stuff
% make sure it's just one channel, find channel

if(iscell(cfg.channel))    % if channel is specified as a cell: need to check that is just one channel
    if(numel(cfg.channel)~=1)
        error('please select one channel to be plotted');
    else
        cfg.channel = cfg.channel{1};
    end
end

channum = find(strcmpi(cfg.channel,data_tf.label)); % only find exact match

dat = data_tf.(cfg.funparameter);

% deal with different dimord setups: perpare data after channel is selected.
chandim = find(ismember(datadimord, 'chan'));

if(chandim==1)   % is there a more elegant way to do this? is it even necessary cause data come from hilbertstats script most likely anyway?!
    dat = squeeze(dat(channum,:,:));
elseif(chandim ==2)
    dat = squeeze(dat(:,channum,:));
else
    dat = squeeze(dat(:,:,channum));
end

%% mess with time:
data_tf.time = data_tf.time*1000;



%% Mask data

% first some actions on the maskparameter field, that are also valied with marking
if(isfield(cfg, 'maskparameter'))
    % if there is a maskparameter field: check its validity
    
    mask = data_tf.(cfg.maskparameter);
    
    % check if maskparameter has right dimensions
    dims = size(data_tf.(cfg.funparameter));
    if(sum(dims~=size(mask))~=0)   % whatever
        error('Maskparameter dimensions do not fit data dimensions.')
    end
    
    
    if(chandim==1)   % is there a more elegant way to do this? is it even necessary cause data come from hilbertstats script most likely anyway?!
        mask = squeeze(mask(channum,:,:));
    elseif(chandim ==2)
        mask = squeeze(mask(:,channum,:));
    else
        mask = squeeze(mask(:,:,channum));
    end
    
    % check if mask is truly binary
    bincheck = unique(mask);
    if(~islogical(bincheck))
        error('Maskparameter is not binary');
    end
    
end

% check if masking and marking is required and throw error
if(cfg.masking && cfg.marking)
    error('Masking and marking at the same time does not make any sense')
end


if(cfg.masking)
    
    if(~isfield(cfg, 'maskparameter'))
        error('No masking parameter specified')
    end
    
    dat(~mask) = NaN;
    
end

%% Check the frequency bands

if(~isfield(data_tf, 'freqbands'))              % if no data_tf.freqbands
    if(isfield(cfg, 'freqbands'))               % but cfg.freqbands
        data_tf.freqbands = cfg.freqbands;      % then add this to data_tf
    else
        error('Please specify cfg.freqbands');  % if neither: quit.
    end
end


% check if freq spec. has 2 dimensions
if(size(data_tf.freqbands,2)~=2)
    error('You have to specify frequency bands as a n_freq x 2 matrix');
end

%% Data dimensions

% dat is now either freq x time or time x freq, we want to have it freq x time
if(size(dat,1)~=size(data_tf.freqbands,1))
    dat = dat';
end


%% finding "clusters" to mark if marking == true

if(cfg.marking)
    
    if(~isfield(cfg, 'maskparameter'))
        error('No masking parameter specified')
    end
    
    % prealloc
    clusterson = cell(1, size(dat,1)); clustersoff = cell(1, size(dat,1));
    
    for ii=1:size(dat,1)
        % this needs to be done for every frequency separately
        maskcumlsum = cumsum(mask(ii,:));
        cluststart = []; countstart = []; cluststop = []; countstop = [];
        
        if(sum(maskcumlsum)==0)   % no significant timepoints at all
            clusterson{ii} =   [];
            clustersoff{ii} =  [];
        else
            
            % now find the "edges" of "clusters"
            countstart = 0; countstop = 0; % counter for start and stop of clusters
            
            % STARTPOINTS OF CLUSTERS
            % first check the 1st sample:
            if(maskcumlsum(1)==1)   % if the cluster starts at the beginning of the whole thing
                % check if cluster is only 1 point
                if(maskcumlsum(2)==1)
                    % do nothing
                else
                    countstart = countstart + 1;
                    cluststart(countstart) = 1;
                end
                
            end
            
            if(maskcumlsum(2)==1 && maskcumlsum(1)==0)
                    countstart = countstart + 1;  
                    cluststart(countstart) = 1;
            end
                
            
            for n=2:length(maskcumlsum)-2 % we do not want to check the last point here
                
                if(maskcumlsum(n) < maskcumlsum(n+1) && maskcumlsum(n) == maskcumlsum(n-1))
                    countstart = countstart + 1;
                    cluststart(countstart) = n+1;
                end
            end
            
            % ENDPOINTS OF CLUSTERS
            for n=2:length(maskcumlsum)-1
                
                if(maskcumlsum(n) > maskcumlsum(n-1) && maskcumlsum(n) == maskcumlsum(n+1))
                    countstop = countstop + 1;
                    cluststop(countstop) = n;
                end
                
            end
            
            % check the last sample:
            if(mask(ii,end)==1)
                % check if cluster is only 1 point
                if(mask(ii,end-1)==0)
                    % do nothing
                else
                countstop = countstop + 1;
                cluststop(countstop) = length(maskcumlsum);
                end
            end
            
            % now check if clusters are sufficently large
            
                
            if(~isfield(cfg, 'markingclust'))
                error('You have to specify how many samples form a cluster in cfg.markingclust.')
            end
            
            if(~isnumeric(cfg.markingclust))
                error('cfg.markingclust should be numeric')
            end
            
            if(cfg.markingclust == 1)
                error('clusters of 1 are not supported')
            end
            
            clustdiff = cluststop-cluststart;
            
            idx = find(clustdiff<cfg.markingclust);
            cluststop(idx)=[]; cluststart(idx)=[];
            
            % samples should be converted to time!
            
            
            clusterson{ii} =   data_tf.time(cluststart);
            clustersoff{ii} =  data_tf.time(cluststop);
            
        end
        % if we do that, we need a copy of the original freqbands (before we maybe mess with it due to ommitted freqs):
        freqorig = data_tf.freqbands;
        
    end
end


%% find possible gaps in frequencybands

for ii=2:size(data_tf.freqbands,1)
    differ(ii-1) = data_tf.freqbands(ii,1)-data_tf.freqbands(ii-1,2);
end

idx = find(differ~=0);

if(size(idx,2)~=0)      % if there are gaps in data
    idx = fliplr(idx);  % flip the index to avoid shifting of rows
    for ii=1:size(idx,2)
        % fill the gap in the freqbands definition
        gapfreq = [data_tf.freqbands(idx(ii),2), data_tf.freqbands(idx(ii)+1,1)];
        data_tf.freqbands = [data_tf.freqbands(1:idx(ii),:); gapfreq; data_tf.freqbands(idx(ii)+1:end,:)];
        
        % fill NaNs in data at appropriate position
        timedim = size(dat,2);
        nandummy = NaN(1,timedim);
        dat = [dat(1:idx(ii),:); nandummy; dat(idx(ii)+1:end,:)];
    end
end

%% lower and upper limits in dat for plotting

for jj=1:size(dat,1)
    datplot(jj*2-1,:) = dat(jj,:);
    datplot(jj*2,:)   = dat(jj,:);
end


%% plot

% take care of the zlims
if(~isfield(cfg, 'zlim'))
    cfg.zlim = 'maxabs'; % use maxabs as default
end

if(isnumeric(cfg.zlim))   % if user specified own values
    climlo = cfg.zlim(1);
    climup = cfg.zlim(2);
else
    
    % make sure zlim is only based on the area that will be plotted
    if(isfield(cfg, 'xlim'))
        xmin = dsearchn(data_tf.time',cfg.xlim(1)); xmax = dsearchn(data_tf.time',cfg.xlim(2));
    else
        xmin = 1; xmax = size(datplot,2);
    end
    
    switch cfg.zlim
        case 'maxabs'
            cmin = min(min(datplot(:,xmin:xmax)));
            cmax = max(max(datplot(:,xmin:xmax)));
            clims = max(abs([cmin, cmax]));    % automatically center zero on color bar
        case 'mincenter'
            clims = abs(min(min(datplot(:,xmin:xmax))));    % centered, max values are min and +min
        case 'maxcenter'
            clims = abs(max(max(datplot(:,xmin:xmax)))); % centered, max values are -max and max
    end
    climlo = -clims;
    climup = clims;
end

% get the freqs in the right order
freqs = reshape(data_tf.freqbands',1,size(data_tf.freqbands,1)*2);

% plot the essential thing:
h=surf(data_tf.time, freqs, datplot);
shading interp


if(cfg.zeroline)
    hold on
        z = round(max(max(datplot))) + 5;   % just make sure that it is "on top" of all the rest;
        plot3([0,0], [freqorig(1,1), freqorig(end,end)], [z,z],...
             'color', [0,0,0], 'linewidth', 1.5); hold on
end
        

% plot marking boxes stuff if needed:
if(cfg.marking)
    z = round(max(max(datplot))) + 5;   % just make sure that it is "on top" of all the rest;
    hold on
    for ii=1:size(mask,1)
        for jj=1:size(clusterson{ii},2)
            plot3([clusterson{ii}(jj), clustersoff{ii}(jj), clustersoff{ii}(jj), clusterson{ii}(jj), clusterson{ii}(jj)], ...   % x axis
                [freqorig(ii,1), freqorig(ii,1),  freqorig(ii,2), freqorig(ii,2), freqorig(ii,1)], ...                          % y axis
                [z, z, z, z, z], ...                                                                                            % z axis
                'color', [0,0,0], 'linewidth', 1.5); hold on
        end
    end
    
    
end

% take care of the rest
view(gca, 2)
title(cfg.channel, 'FontSize', 18);
ylabel('Frequency (Hz)');
xlabel('Time (ms)');
set(gca, 'clim', [climlo climup]);

if(isfield(cfg, 'xlim'))
    set(gca, 'xlim', [cfg.xlim(1), cfg.xlim(2)]);
else
    axis tight
end
grid off
colorbar;

