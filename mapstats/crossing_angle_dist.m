
function [angles, counts, crossing_angles] = crossing_angle_dist(map1, map1_type, map2, map2_type, mask, nangles, bins)
%CROSSING_ANGLE_DIST Calculate contour crossing angle distribution.
%
% [angles, counts, crossing_angles] = crossing_angle_dist(map1, map1_type, map2, map2_type, mask, nangles, bins)
%
% Calculates a distribution of intersection angles between the given two
% maps, which can be any combination of OP, OD, and SF maps. The *_type
% inputs should be set to two character strings, one of 'op', 'od', or 'sf'
% as appropriate. Uses one of od_op_crossing, od_sf_crossing, or
% sf_op_crossing to find the crossing angles. A binary mask can be used
% to restrict the analysis. Both the distribution, raw counts in each bin,
% and all the crossing angles are returned. [0:10:80, 91] is the default
% bins, where the last element is used to correctly handle angles wrapping
% around.
%
% See also:
% od_op_crossing
% od_sf_crossing
% sf_op_crossing

% Handle input
if nargin < 4
    error('At least four inputs are required.');
end
if nargin < 5
    mask = true(size(op));
end
if nargin < 6
    nangles = 4;
end
if nargin < 7
    bins = [0:10:80, 91];
end

% Mask
map1(~mask) = nan;
map2(~mask) = nan;

if strcmp(map1_type, 'op')
    if strcmp(map2_type, 'od')
        % Calculate crossing angles of OP and OD
        crossing_angles = od_op_crossing(map2, map1, nangles);
    elseif strcmp(map2_type, 'sf')
        % Calculate crossing angles of OP and SF
        crossing_angles = sf_op_crossing(map2, map1, nangles);
    end
elseif strcmp(map1_type, 'od')
    if strcmp(map2_type, 'op')
        % Calculate crossing angles of OD and OP
        crossing_angles = od_op_crossing(map1, map2, nangles);
    elseif strcmp(map2_type, 'sf')
        % Calculate crossing angles of OD and SF
        crossing_angles = od_sf_crossing(map1, map2);
    end
elseif strcmp(map1_type, 'sf')
    if strcmp(map2_type, 'op')
        % Calculate crossing angles of SF and OP
        crossing_angles = sf_op_crossing(map1, map2, nangles);
    elseif strcmp(map2_type, 'od')
        % Calculate crossing angles of SF and OD
        crossing_angles = od_sf_crossing(map2, map1);
    end
end

% Calculate distribution
a = histc(crossing_angles*180/pi, bins);
n = a(1:end-1);
counts = n;
angles = n/sum(n);
