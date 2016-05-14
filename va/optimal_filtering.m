
function [f_od, f_op, a] = optimal_filtering(t_od, t_op, od, op)
%OPTIMAL_FILTERING Low- and high-pass the given OP and OD maps optimally.
%
% [f_od, f_op, a] = optimal_filtering(t_od, t_op, od, op)
%
% Low- and high-pass Gaussian filters the maps od and op so as to maximise
% their linear correlation with the maps t_od and t_op. The filtered maps,
% as well as the wavelengths of the two optimal filters a, are returned.

fun = @(a)(1 - filtered_correlation(a, t_od, t_op, od, op));
a0 = [5, 10];
options = optimset('Display', 'None');
a = fminsearch(fun, a0, options);

[f_od, f_op] = filter_maps(a, od, op);

end


function [f_od, f_op] = filter_maps(a, od, op)
% Gaussian filter the maps od and op with wavelengths given by a.

a1size = max([10, ceil(6*a(1))]);
a2size = max([10, ceil(6*a(2))]);
if a(1) > 0
    a1filter = fspecial('Gaussian', a1size, a(1));
end
if a(2) > 0
    a2filter = fspecial('Gaussian', a2size, a(2));
end

if a(2) > 0
    f_od = od - imfilter(od, a2filter, 'symmetric');
else
    f_od = od;
end
if a(1) > 0
    f_od = imfilter(f_od, a1filter, 'symmetric');
end

if a(2) > 0
    f_op = op - imfilter(op, a2filter, 'symmetric');
else
    f_op = op;
end
if a(1) > 0
    f_op = imfilter(f_op, a1filter, 'symmetric');
end

end


function cc = filtered_correlation(a, t_od, t_op, od, op)
% Calculate the linear correlation between t_od, t_op and od, op after
% od and op are filtered by the parameters given in a.

[f_od, f_op] = filter_maps(a, od, op);

filtered_maps = vertcat(f_od(:), real(f_op(:)), imag(f_op(:)));
true_maps = vertcat(t_od(:), real(t_op(:)), imag(t_op(:)));

cc = corr(filtered_maps, true_maps);

end
