function mat = nemo_ftsource2mat(source,subfield)
inside_idx = find(source.inside);
if(isfield(source,subfield))
    mat = reshape([source.(subfield){:}],length(source.freq),length(source.time),length(inside_idx));
elseif(isfield(source.avg,subfield))
    mat = reshape([source.avg.(subfield){:}],length(source.freq),length(source.time),length(inside_idx));
else
    error([subfield ' not found in structure']);
end
mat = permute(mat, [3 1 2]);