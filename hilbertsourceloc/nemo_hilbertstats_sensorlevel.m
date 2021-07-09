%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% perform stats on Hilbert sensorlevel data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should best be used with nemo_hilbert_sensorlevel


baselinewindow_idx = dsearchn(data_hilb_bp{ii}.time{1}',baselinewindow');    % this is a hack that only works if all trials have same time axis! FIX THIS.

%%
disp(freqbands(ii,:));
aph = zeros(length(data_hilb_bp{ii}.trial), numel(data_hilb_bp{ii}.label), length(data_hilb_bp{ii}.time{1}));
aa = aph;
for jj=1:length(data_hilb_bp{ii}.trial)
    trials(jj,:,:) = data_hilb_bp{ii}.trial{jj};
    aph(jj,:,:)    = angle(data_hilb_bp{ii}.trial{jj});
    aa(jj,:,:)     = abs(data_hilb_bp{ii}.trial{jj});
end


phvar=squeeze(circ_var(aph,[],[],1));

aamean=squeeze(mean(aa,1)); % Hilbert analytic amplitude
aabaseline = mean(aamean(:,baselinewindow_idx(1):baselinewindow_idx(2)),2);
aameannorm = 20*log10(aabaseline./aamean);  % normalize against baseline

% create a dummy tlk structure to use for output storing
data_hilbert{ii} = ft_timelockanalysis([], data_hilb_bp{ii}); % initialize source_hilb
data_hilbert{ii}.avg = aameannorm;
data_hilbert{ii}.itc = data_hilbert{ii}.avg;
data_hilbert{ii}.itc = 1 - phvar;


if(tfstats)
    % TODO: FDR correction? cluster analysis?
    p = ones(size(aa,2),size(aa,3));
    zval = zeros(size(aa,2),size(aa,3));
    %     pks = p; ksval = zval;
    
    ft_progress('init','etf');
    for jj=1:size(aa,2)
        aabl = squeeze(aa(:,jj,baselinewindow_idx(1):baselinewindow_idx(2)));
        aabl = aabl(:);
        for kk=1:size(aa,3)
            [p(jj,kk),~,stats]=ranksum(aa(:,jj,kk),aabl);
            zval(jj,kk)=stats.zval;
            %             [~,pks(jj,kk),ksval(jj,kk)] = kstest2(aa(:,jj,kk),aabl);
        end
        ft_progress(jj/size(aa,2),'%d of %d',jj,size(aa,2));
    end
    ft_progress('close');
    
    data_hilbert{ii}.stat = zval;
    data_hilbert{ii}.pval = log10(p);
    
end
