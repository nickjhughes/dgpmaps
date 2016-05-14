
function crossing = od_sf_crossing(od, sf)
%OD_SF_CROSSING Calculate angles of intersection of SF and OD contours.
%
% crossing = od_sf_crossing(od, sf)
%
% Calculates the angle of intersection of all crossings of the iso-contours of
% the given SF and OD maps. The zero-level contour of the OD map is used,
% and 4 log-spaced contours between log10(0.5) and log10(2) are used for
% the SF map.
%
% Requires intersections:
% http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections/content/intersections.m
%
% See also:
% crossing_angle_dist

% Calculate contours
od_c = contourc(od, [0 0]);
sf_vals = logspace(log10(0.5), log10(2), 4);
sf_c = contourc(sf, sf_vals);

% Correct contours
od_c(:,od_c(1,:)==0) = nan;
sf_c(:,sf_c(1,:)==0) = nan;

% Initialise contours
x1 = od_c(1,:);
y1 = od_c(2,:);
x2 = sf_c(1,:);
y2 = sf_c(2,:);

% Make sure there are contours
if isempty(x1) || isempty(x2)
    crossing = [];
    return
end

% Find intersections
[x0, ~, I, J] = intersections(x1, y1, x2, y2);

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
