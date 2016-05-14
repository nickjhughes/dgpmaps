
function plot_filter(filter, filename, cmap_filter, figure_width, figure_height)
%PLOT_FILTER Plot the given filtering kernel in 3D.
%
% plot_filter(filter, filename, cmap_filter, figure_width, figure_height)
%
% Plots the given filtering kernel in 3D, saves the image if a filename
% is given, using the given colormap (default is HSV), with figure widths
% and heights given in centimeters (defaults are 8 and 6).
%
% Requires export_fig:
% http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

if nargin < 2
    filename = nan;
end
if nargin < 3
    cmap_filter = hsv(16);
end
if nargin < 5
    figure_width = 8;
    figure_height = (3/4)*figure_width;
end

figure;
surf(filter, 'LineStyle', 'None');
view(3);
colormap(cmap_filter);
axis off;
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0.1, 0.1, -0.1, -0.1]);
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Centimeters', 'Position', [10 10 figure_width figure_height]);
if ~isnan(filename)
    export_fig(filename);
    close;
end
