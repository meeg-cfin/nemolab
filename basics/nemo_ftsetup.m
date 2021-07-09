codepath = '~/code'; % location of your Matlab toolboxes
addpath([codepath '/fieldtrip'])
%% set up FieldTrip (add paths etc.)
ft_defaults
ftpath = fileparts(which('ft_defaults'));
addpath(genpath([ftpath '/template/atlas/']));
%rmpath(genpath([ftpath '/external']));
addpath([ftpath '/contrib/nutmegtrip/']);

%% add extra toolboxes
addpath([codepath '/nemolab/basics'])
addpath([codepath '/cmocean'])
addpath([codepath '/spm12'])
addpath([codepath '/circstat-matlab/']);



