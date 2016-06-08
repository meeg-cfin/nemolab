% perform stats on Hilbert source data
baselinewindow_idx = dsearchn(source_hilbtmp{1}.time',baselinewindow');
inside_idx = find(source_hilbtmp{1}.inside);

%%
disp(freqbands(ii,:));
aph = zeros(length(source_hilbtmp{ii}.trial),length(inside_idx),length(source_hilbtmp{ii}.trial(1).mom{inside_idx(1)}));
aa = aph;
for jj=1:length(source_hilbtmp{ii}.trial)
    trials(jj,:,:) = cell2mat(source_hilbtmp{ii}.trial(jj).mom);
    aph(jj,:,:) = angle(cell2mat(source_hilbtmp{ii}.trial(jj).mom));
    aa(jj,:,:) = abs(cell2mat(source_hilbtmp{ii}.trial(jj).mom));
end
if(saveRAM)
    source_hilbtmp{ii} = []; % don't use "clear" since this is a for loop!!
end

phvar=squeeze(circ_var(aph,[],[],1));

aamean=squeeze(mean(aa,1)); % Hilbert analytic amplitude
aabaseline = mean(aamean(:,baselinewindow_idx(1):baselinewindow_idx(2)),2);
aameannorm = 20*log10(diag(1./aabaseline)*aamean);  % normalize against baseline

source_hilb{ii} = source_bp{ii}; % initialize source_hilb
source_hilb{ii}.avg.itc = source_hilb{ii}.avg.mom; % initialize itc
for jj=1:length(inside_idx)
    source_hilb{ii}.avg.itc{inside_idx(jj)} = 1 - phvar(jj,:);
    source_hilb{ii}.avg.mom{inside_idx(jj)} = aameannorm(jj,:);
end

if(tfstats)
    % TODO: FDR correction? cluster analysis?
    p = ones(size(aa,2),size(aa,3));
    zval = zeros(size(aa,2),size(aa,3));
%     pks = p; ksval = zval;
    
    ft_progress('init','etf');
    for jj=1:size(aa,2)
        aabl = squeeze(aa(:,jj,baselinewindow_idx(1):baselinewindow_idx(2)));
        aabl = aabl(:);
        parfor kk=1:size(aa,3)
            [p(jj,kk),~,stats]=ranksum(aa(:,jj,kk),aabl);
            zval(jj,kk)=stats.zval;
%             [~,pks(jj,kk),ksval(jj,kk)] = kstest2(aa(:,jj,kk),aabl);
        end
        ft_progress(jj/size(aa,2),'%d of %d',jj,size(aa,2));
    end
    ft_progress('close');
    
    source_hilb{ii}.stat = source_hilb{ii}.avg.mom; % initialize with 'mom' since it's the same size
    source_hilb{ii}.pval = source_hilb{ii}.stat;
    
    for jj=1:length(inside_idx)
        source_hilb{ii}.stat{inside_idx(jj)} = zval(jj,:);
        source_hilb{ii}.pval{inside_idx(jj)} = log10(p(jj,:));
    end
end
