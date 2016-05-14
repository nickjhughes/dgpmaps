
function [pinwod_dist, counts] = pinwod(op, od, mask, nbins, rpinw)
%PINWOD Calculate histogram of pinwheel/OD spatial relationship.
%
% [pinwod_dist, counts] = pinwod(op, od, mask, nbins, rpinw)
%
% Calculates a distribution of OP pinwheel locations relative to an OD map.
% The positive and negative parts of the OD map are split into nbins (5 by
% default) bins such that each contains the same number of pixels. Each OP
% pinwheel is assigned to the appropriate bin, and the positive and negative
% bins are combined. A binary mask can be provided to restrict the analysis.
% If repeatedly analysing a pair of maps, pinwheel location can be pre-computed
% and provided as rpinw for efficiency. The second output contains the raw
% pinwheel counts in each bin.
%
% This metric is a modification of the metric used in:
% Hubener M, Shoham D, Grinvald A, Bonhoeffer T (1997) "Spatial relationships
% among three columnar systems." J Neurosci 17(23):9270-9284.

% Input defaults and validation
if nargin < 2
    error('At least two inputs are required.');
end
if nargin < 3
    mask = true(size(op));
end
if nargin < 4
    nbins = 5;
end
if nargin < 5
    rpinw = [];
end

% Mask
op(~mask) = nan;
od(~mask) = nan;

% Split the OD map into positive and negative sections
od_pos = od;
od_pos(od <= 0) = nan;
od_neg = od;
od_neg(od > 0) = nan;

% Divide the two OD map parts into nbins quantiles
% Positive
if nbins == 2
    q = quantile(od_pos(:), 2);
    q = q(1);
else
    q = quantile(od_pos(:), nbins-1);
end
bin_od_pos = nan(size(od_pos));
bin_od_pos(od_pos < q(1)) = nbins;
for k = 2:nbins-1
    bin_od_pos(od_pos < q(k) & od_pos > q(k-1)) = nbins - k + 1;
end
bin_od_pos(od_pos > q(nbins-1)) = 1;
% Negative
if nbins == 2
    q = quantile(od_neg(:), 2);
    q = q(1);
else
    q = quantile(od_neg(:), nbins-1);
end
bin_od_neg = nan(size(od_neg));
bin_od_neg(od_neg < q(1)) = nbins;
for k = 2:nbins-1
    bin_od_neg(od_neg < q(k) & od_neg > q(k-1)) = nbins - k + 1;
end
bin_od_neg(od_neg > q(nbins-1)) = 1;

% Find pinwheels in OP
if isempty(rpinw)
    rpinw = locate_pinwheels(op);
end
npinw = size(rpinw, 1);

% Calculate histogram of pinwheels in OD quantiles
pinwod_dist = zeros(2, nbins);
for l = 1:npinw
    % Determine if its in the negative or positive section
    if isnan(bin_od_pos(round(rpinw(l,2)), round(rpinw(l,1))))
        % Negative
        p = bin_od_neg(round(rpinw(l,2)), round(rpinw(l,1)));
        if ~isnan(p)
            pinwod_dist(2,p) = pinwod_dist(2,p) + 1;
        end
    else
        % Positive
        p = bin_od_pos(round(rpinw(l,2)), round(rpinw(l,1)));
        if ~isnan(p)
            pinwod_dist(1,p) = pinwod_dist(1,p) + 1;
        end
    end
end

% Combine the two histograms
pinwod_dist = pinwod_dist(1,:) + fliplr(pinwod_dist(2,:));
counts = pinwod_dist;
pinwod_dist = pinwod_dist/sum(pinwod_dist);
