% set up FieldTrip (add paths etc.)

ft_defaults
ftpath = fileparts(which('ft_defaults'));
addpath(genpath([ftpath '/template/atlas/']));
rmpath(genpath([ftpath '/external']));
addpath([ftpath '/contrib/nutmegtrip/']);

%% specific to Sarang's setup
addpath('/Data/bem/spm8newseg');
switch(computer)
    case 'PCWIN64'
        addpath('c:\Users\sarang\Desktop\lab\MATLAB\CircStat2012a');
    case 'GLNXA64'
        addpath(genpath('/Data/MATLAB/stable/CircStat2012a/'));
end
