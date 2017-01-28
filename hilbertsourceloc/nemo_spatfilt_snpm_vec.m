function [stats,w] = nemo_spatfilt_snpm(cfg,dat,eta);
% cfg requires:
%      method = 'lcmv' or 'minnorm'
%      stat = 'ranksum' 'ranksumactive' 'signrank' 'signrankactive'
% dat requires:
%      b
%      L
%      invCy
%      cfg.t_idx (for ranksum or signrank)
%      controlwin
%      activewin


% blsel = 1:201;
% actsel = 231:431;

lf = dat.L * eta;

switch(cfg.method)
    case 'lcmv'
        w = inv(lf' * dat.invCy * lf) * lf' * dat.invCy;

    case 'lcmvvec'
        W = inv(dat.L' * dat.invCy * dat.L) * dat.L' * dat.invCy;
    case 'lcmvvec2'
        for ii=1:3
            W2(ii,:) = inv(dat.L(:,ii)' * dat.invCy * dat.L(:,ii)) * dat.L(:,ii)' * dat.invCy;
        end
    case 'minnorm'
        G = lf * lf'; % Gram matrix
        lambda = min(eig(dat.Cy));
        invG = inv(G + lambda * eye(size(G))); % regularized G^-1
        w = (invG * lf)';
end
 

switch(cfg.stat)
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
    case 'ranksum'
        aa = abs(s(cfg.t_idx,:));
        [pval,~,stats]=ranksum(aa,aabl(:));
    case 'ranksumactive'
        aabl = abs(w * dat.bcontrol);
        aaact = abs(w * dat.bactive);
        [pval,~,stats]=ranksum(aaact,aabl);
%        [~,~,~,stats]=ttest2(aaact(:),aabl(:));
    case 'power'
%        Sact = mean((w * real(dat.bactive)).^2);
%        Scon = mean((w * real(dat.bcontrol)).^2);

        Sact = w * dat.Ract * w';
        Scon = w * dat.Rcon * w';
        pval = 0;
        stats.zval = log10(Sact/Scon);
%        [~,~,~,stats]=ttest2(aaact(:),aabl(:));
    case 'median'
      Sact = median((w * dat.bact).^2);
      Scon = median((w * dat.bcon).^2);
 

        pval = 0;
        stats.zval = log10(Sact/Scon);
%        [~,~,~,stats]=ttest2(aaact(:),aabl(:));
end


stats.pval = pval;
