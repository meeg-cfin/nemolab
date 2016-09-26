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
%   cfg.maskparameter      = if masking is desired, this should point to a field in data_tf that contains a binary mask
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




%% prepare data
% find channel
channum = strmatch(cfg.channel,data_tf.label);

dat = data_tf.(cfg.funparameter);

% check which dimension is channel
datadimord = tokenize(data_tf.dimord, '_');

if(length(datadimord)~=3)
    error('can only handle 3-dimensional data')
end

chandim = find(ismember(datadimord, 'chan'));

% deal with different dimord setups
if(chandim==1)   % is there a more elegant way to do this? is it even necessary cause data come from hilbertstats script most likely anyway?!
    dat = squeeze(dat(channum,:,:));
elseif(chandim ==2)
    dat = squeeze(dat(:,channum,:));
else
    dat = squeeze(dat(:,:,channum));
end


%% Mask data

if(isfield(cfg, 'maskparameter'))
    
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

%% find possible gaps in frequencybands
for ii=2:size(data_tf.freqbands,1)
   differ(ii-1) = data_tf.freqbands(ii,1)-data_tf.freqbands(ii-1,2);
end

idx = find(differ~=0);

if(size(idx,2)~=0)   % if there are gaps in data
    idx = fliplr(idx);  % flip the index to avoid shifting rows
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

if(~isfield(cfg, 'zlim'))
    cfg.zlim = 'maxabs'; % use maxabs as default
end

if(isnumeric(cfg.zlim))   % if user specified own values
    climlo = cfg.zlim(1);
    climup = cfg.zlim(2);
else  
    switch cfg.zlim
        case 'maxabs'
            cmin = min(min(datplot));
            cmax = max(max(datplot));
            clims = max(abs([cmin, cmax]));    % automatically center zero on color bar
        case 'mincenter'
            clims = abs(min(min(datplot)));    % centered, max values are min and +min
        case 'maxcenter'
            clims = abs(max(max(datplot))); % centered, max values are -max and max
    end
    climlo = -clims;
    climup = clims;
end

freqs = reshape(data_tf.freqbands',1,size(data_tf.freqbands,1)*2);

h=surf(data_tf.time, freqs, datplot);
shading interp
view(gca, 2)
title(cfg.channel, 'FontSize', 18);
ylabel('Frequency (Hz)');
xlabel('Time (s)');
set(gca, 'clim', [climlo climup]);
if(isfield(cfg, 'xlim'))
    set(gca, 'xlim', [cfg.xlim(1), cfg.xlim(2)]);
end
axis tight
grid off
colorbar;
