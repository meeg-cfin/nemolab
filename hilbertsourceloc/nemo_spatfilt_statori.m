function [ori] = nemo_spatfilt_statori(cfg,dat);
%function [ori,s,stat] = nemo_spatfilt_statori(cfg,dat);

orires = pi/32;
[th,ph]=meshgrid(0:orires:pi,(-pi/2):orires:(pi/2));
[etax,etay,etaz] = sph2cart(th,ph,1);
eta = [etax(:) etay(:) etaz(:)]';

%%
clear stat
%dat.W = inv(dat.L' * dat.invCy * dat.L) * dat.L' * dat.invCy;
%parfor ii=1:size(eta,2)
for ii=1:size(eta,2)
     [stat(ii)] = nemo_spatfilt_snpm(cfg,dat,eta(:,ii));
end
zval = [stat(:).zval];
%%
%zval = [stat(:).tstat];
[maxzval,phsel] = max(abs(zval));
ori = eta(:,phsel); %/norm(eta(:,phsel));
%zval = zval(phsel)

%[s,stat] = nemo_spatfilt_snpm(cfg,dat,ori);