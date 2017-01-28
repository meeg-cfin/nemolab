%%
nemo_ftsetup % set up path etc.

basepath = pwd;

mripath = [basepath '/s' cfgnemo.participant '.nii'];
normmripath = [basepath '/ws' cfgnemo.participant '.nii'];
% mri = ft_read_mri(mripath,'dataformat','nifti');
% mri.coordsys = 'spm';

cfg = [];
cfg.atlas = ft_read_atlas('ROI_MNI_V4.nii'); % AAL atlas
cfg.mripath = normmripath;

inside_idx = find(source_tf.inside);

if(~exist('cfgnemo','var'))
    cfgnemo.sourceplottype = 'fancy';
end

switch(cfgnemo.sourceplottype)
    case 'simple'
        %% simple plot
        cfg.mripath = normmripath;
        cfg.funparameter = 'avg.mom';
        cfg.evokedoverlay = 1;
        % cfg.plottype = 'ts';  % if you want to view as time series instead
        nmt_sourceplot(cfg,source_tf);
        
    case 'fancy'
        %% fancy plotting with thresholding (e.g. with p-value)
        cfg.funparameter = {'itc','stat'};
        cfg.maskparameter = 'msk';
        masktype = 'pval';  % masktypes: pval, momhigh, momlow, momabs, itchigh, itclow, pow
        maskthresh = log10(.001); for pval masktype, this should be log(pval)
        
        if(isfield(source_tf,'dim'))
            % FIXME: source_tf_avg contains "dim", but others don't. this confuses
            % source_tfplot... long-term solution needed, but this is a temporary fix:
            source_tf = rmfield(source_tf,'dim');
        end
        
        if(0)
            switch(cfg.funparameter) % works for pow but not mom...
                case 'avg.pow'
                    %        cfg.maskparameter = 'msk';
                case 'avg.itc'
                    itc = reshape([source_tf.avg.itc{:}],length(source_tf.time),length(inside_idx));
                    
                    
                    
                case 'avg.mom' % FIXME: this is really messy, because of the cell strucutre of avg.mom :-(
                    % nmt_polaritytweak;
                    mom = reshape([source_tf.avg.mom{:}],length(source_tf.time),length(inside_idx));
                    [peakval,peakvox]=max(abs(mom),[],2);
                    
                case 'avg.mompow'
                    source_tf.msk = zeros(size(source_tf.avg.mom,1),1,length(source_tf.time));
                    for tt=1:length(source_tf.time) % pick threshold relative to peak at each given time point
                        source_tf.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*max(abs(mom(:)));
                    end
            end
        end
        
        source_tf.msk = zeros(size(source_tf.avg.mom,1),length(source_tf.freq),length(source_tf.time));
        switch(masktype)
            case 'itchigh'
                itc = reshape([source_tf.avg.itc{:}],length(source_tf.freq),length(source_tf.time),length(inside_idx));
                itc = permute(itc, [3 1 2]);
                source_tf.msk(inside_idx,1:length(source_tf.freq),:) = itc > maskthresh;
            case 'itclow'
                itc = reshape([source_tf.avg.itc{:}],length(source_tf.freq),length(source_tf.time),length(inside_idx));
                itc = permute(itc, [3 1 2]);
                source_tf.msk(inside_idx,1:length(source_tf.freq),:) = itc < maskthresh;
            case 'pow'
                source_tf.msk = source_tf.avg.pow > maskthresh*max(source_tf.avg.pow); % threshold at some % of peak
            case 'momhigh'
                source_tf.msk(inside_idx,1,:) = mom' > maskthresh;
            case 'momlow'
                source_tf.msk(inside_idx,1,:) = mom' < maskthresh;
            case 'momabs'
                source_tf.msk(inside_idx,1,:) = abs(mom') > maskthresh;
            case 'momratio'
                for tt=1:length(source_tf.time) % pick threshold relative to peak at each given time point
                    %            source_tf.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*peakval(tt);
                    source_tf.msk(inside_idx,1,tt) = abs(mom(tt,:))>maskthresh*max(abs(mom(:)));
                end
            case 'pval'
                pval = reshape([source_tf.pval{:}],length(source_tf.freq),length(source_tf.time),length(inside_idx));
                pval = permute(pval, [3 1 2]);
                source_tf.msk(inside_idx,1:length(source_tf.freq),:) = pval < maskthresh;
            otherwise
                error('masktype unknown');
        end
        
        
        nmt_sourceplot(cfg,source_tf);
        
        
end

