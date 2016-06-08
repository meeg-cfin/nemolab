function source_out = ft_apply_spatialfilter(dat,source_in)
% TODO
% - allow for 3 orientations when supplied with averaged data
% - allow for 2 orientations for oldschool fun
% - wrap 3-orientation output in fieldtrip's expected style

source_out = source_in; % initialize source_out
inside_idx = find(source_in.inside)';
spatfilt = vertcat(source_in.avg.filter{inside_idx});


if(isfield(dat,'avg'))
    mom = spatfilt * dat.avg;
    
    for i=1:length(inside_idx)
        source_out.avg.mom{inside_idx(i)} = mom(i,:);
    end
elseif(isfield(dat,'trial'))   % calculate the virtual sensors for more than one trial
    num_oris = size(source_in.avg.filter{inside_idx(1)},1);
    
    chansel = match_str(dat.label,source_in.avg.label); % select only channels that went into spatial filter construction
    
    if(num_oris == 1)
        spatfilt=cell2mat(source_in.avg.filter);
        
        for ii=1:size(dat.trial,2)
            source_out.trial(ii).mom = cell(size(source_in.avg.mom)); % initialize mom cells
            s = spatfilt*dat.trial{ii}(chansel,:);
            source_out.trial(ii).mom(inside_idx) = mat2cell(s,ones(size(s,1),1),size(s,2));
        end
        
        
        
    elseif(num_oris ==3)
        
        
        spatfilt = reshape([source_in.avg.filter{inside_idx}],3, numel(source_in.avg.label),length(inside_idx) );
        
        for ii=1:length(dat.trial)
            virtsens.x{ii} = squeeze(spatfilt(1,:,:))'*dat.trial{ii}(chansel,:);
            virtsens.y{ii} = squeeze(spatfilt(2,:,:))'*dat.trial{ii}(chansel,:);
            virtsens.z{ii} = squeeze(spatfilt(3,:,:))'*dat.trial{ii}(chansel,:);
        end
        
        source_out.virtsens = virtsens;      % TO DO: could of course be shaped in a nicer way, e.g. in one cell
        
    else
        error('Cannot handle %d orientations yet, sorry.', num_oris)   % TO DO
    end
    
    % not yet implemented
    source_out.avg.pow = [];
    
end
