
function [od, op] = vector_avg(data, eyes, orientations)
%VECTOR_AVG Calculate vector averaged OP and OD maps.
%
% [od, op] = vector_avg(data, eyes, orientations)
%
% Calculates vector averaged OP and OD maps from the given data array,
% which should be of size [height, width, num_eyes, num_orientations, 
% num_trials]. If monocular data is used eyes should be [-1, 1], and
% orientations should contain the stimulus orientations in degrees.

% Inputs
if nargin < 2
    eyes = [-1, 1];
end
if nargin < 3
    orientations = [0, 45, 90, 135];
end

% Vector average
data = mean(data,5);
binoc_sing_conds = squeeze(mean(data,3));
op = zeros(size(data,1), size(data,2));
for j = 1:length(orientations)
    op = op + binoc_sing_conds(:,:,j)*cosd(2*orientations(j)) + 1i*binoc_sing_conds(:,:,j)*sind(2*orientations(j));
end
iso_sing_conds = squeeze(mean(data,4));
od = eyes(2)*iso_sing_conds(:,:,2) + eyes(1)*iso_sing_conds(:,:,1);
