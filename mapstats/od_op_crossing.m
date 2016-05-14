
function crossing = od_op_crossing(od, op, nangles)
%OD_OP_CROSSING Calculate angles of intersection of OP and OD contours.
%
% crossing = od_op_crossing(od, op, nangles)
%
% Calculates the angle of intersection of all crossings of the iso-contours of
% the given OP and OD maps. The zero-level contour of the OD map is used,
% and nangles (default is 4) equally spaced iso-orientation contours of the OP
% map.
%
% Requires intersections:
% http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections/content/intersections.m
%
% See also:
% crossing_angle_dist

% Input defaults
if nargin < 3
    nangles = 4;
end

% Initialise
or_real = abs(op).*cos((nangles/2)*angle(op)/2);
or_imag = abs(op).*sin((nangles/2)*angle(op)/2);

% Calculate contours
od_c = contourc(od, [0 0]);
or_real_c = contourc(or_real, [0 0]);
or_imag_c = contourc(or_imag, [0 0]);

% Combine contours
od_c(:,od_c(1,:)==0) = nan;
or_real_c(:,or_real_c(1,:)==0) = nan;
or_imag_c(:,or_imag_c(1,:)==0) = nan;
or_c = [or_real_c, [nan; nan], or_imag_c];

% Initialise contours
x1 = od_c(1,:);
y1 = od_c(2,:);
x2 = or_c(1,:);
y2 = or_c(2,:);

% Make sure there are contours
if isempty(x1) || isempty(x2)
    crossing = [];
    return
end

% Find intersections
[x0, y0, I, J] = intersections(x1, y1, x2, y2);
[x0_p, y0_p] = intersections(or_real_c(1,:), or_real_c(2,:), or_imag_c(1,:), or_imag_c(2,:));

% Remove intersections that are less than half a pixel away from a pinwheel
remove = nan(size(x0,1),1);
r = 1;
for j = 1:size(x0,1)
    for k = 1:size(x0_p,1)
        if sqrt((x0(j) - x0_p(k))^2 + (y0(j) - y0_p(k))^2) < 0.5
            remove(r) = j;
            r = r + 1;
        end
    end
end
remove = remove(1:r-1);
x0(remove) = [];
y0(remove) = [];
I(remove) = [];
J(remove) = [];

% Calculate crossing angles
crossing = nan(length(x0),1);
for j = 1:length(x0)
    x = [x1(floor(I(j))),x1(floor(I(j))+1)];
    y = [y1(floor(I(j))),y1(floor(I(j))+1)];
    th1 = atan2(y(2)-y(1),x(2)-x(1));
    x = [x2(floor(J(j))),x2(floor(J(j))+1)];
    y = [y2(floor(J(j))),y2(floor(J(j))+1)];
    th2 = atan2(y(2)-y(1),x(2)-x(1));
    angle_diff = abs(th1 - th2);
    angle_diff = min(min(angle_diff, abs(pi-angle_diff)), 2*pi-angle_diff);
    crossing(j) = angle_diff;
end
