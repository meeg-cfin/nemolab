function [stats,w] = nemo_spatfilt_snpm(cfg,dat,eta);
% cfg requires:
%      method = 'lcmv' or 'minnorm'
%      stat = 'ranksum' 'ranksumactive' 'signrank' 'signrankactive'
% dat requires:
%      L
%      invC

%      b (for signrank*)
%      cfg.t_idx (for ranksum or signrank)
%      controlwin (for signrankactive)
%      activewin (for signrankactive)


lf = dat.L * eta;


switch(cfg.method)
    case 'lcmv'
        w = inv(lf' * dat.invC * lf) * lf' * dat.invC;

    case 'lcmvvec'
            W = inv(dat.L' * dat.invC * dat.L) * dat.L' * dat.invC;
            w = eta' * W;
    case 'lcmvvecpre'
            w = eta' * dat.W;
    case 'lcmvvec2'
            W = trace(inv(dat.L' * dat.invC * dat.L)) * dat.L' * dat.invC;
            w = eta' * W;
    case 'lcmvvec3'
        for ii=1:3
            W(ii,:) = inv(dat.L(:,ii)' * dat.invC * dat.L(:,ii)) * dat.L(:,ii)' * dat.invC;
        end
        w = eta' * W;
    case 'minnorm'
        error('NO!!!! THIS IS FUCKED! GRAM MATRIX IS NOT CORRECT!!')
        G = lf * lf'; % Gram matrix
        lambda = min(eig(dat.Cy));
        invG = inv(G + lambda * eye(size(G))); % regularized G^-1
        w = (invG * lf)';
end
 

switch(cfg.stat)
    case 'nai'
            wnai = w / sqrt(w*w');
            stats.zval = wnai * dat.C * wnai';
            pval = 0;
    case 'powrat'
        stats.zval = (w * dat.Cact * w')/(w * dat.Ccon * w');
        pval = 0;
    case 'ranksumactive'
        aabl = abs(w * dat.bcontrol);
        aaact = abs(w * dat.bactive);
        [pval,~,stats]=ranksum(aaact,aabl);
%        [~,~,~,stats]=ttest2(aaact(:),aabl(:));
    case 'power'
%        Sact = mean((w * real(dat.bactive)).^2);
%        Scon = mean((w * real(dat.bcontrol)).^2);

        Sact = w * dat.Cact * w';
        Scon = w * dat.Ccon * w';
        pval = 0;
        stats.zval = log10(Sact/Scon);
%        [~,~,~,stats]=ttest2(aaact(:),aabl(:));
    case 'signrank'
        aa = abs(s(cfg.t_idx,:));
        [pval,~,stats]=signrank(aa,median(aabl));
    case 'signrankactivenew'
        aaact = zeros(1,size(dat.b,3));
        aabl = aaact;
        for jj=1:size(dat.b,3)
            aaact(jj) = w * dat.Ract(:,:,jj) * w';
            aabl(jj) = w * dat.Rcon(:,:,jj) * w';
        end
        [pval,~,stats]=signrank(aaact,aabl);
    case 'signrankactive'
        for jj=1:size(dat.b,3)
            s(:,jj) = w * dat.b(:,:,jj);
        end


        aabl = abs(s(dat.controlwin,:));
        aaact = abs(s(dat.activewin,:));
        [pval,~,stats]=signrank(median(aaact),median(aabl));
    case 'median'
      Sact = median((w * dat.bact).^2);
      Scon = median((w * dat.bcon).^2);
 

        pval = 0;
        stats.zval = log10(Sact/Scon);
    case 'ssd'
        s = abs(w * dat.erm);
        stats.zval =c;
        pval = stats.zval;
end


stats.pval = pval;
