function bw_plot3(pos, args)
% makes nasty repetitions in plot3 function call obsolete

if(~exist('args', 'var'))
    args = {'b.'};
end

plot3(pos(:,1), pos(:,2), pos(:,3), args{:});

axis equal
end
