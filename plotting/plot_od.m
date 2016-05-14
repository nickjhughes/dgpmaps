
function plot_od(od, filename, cmap_od, figure_width, figure_height)
%PLOT_OD Plot the given OD map.
%
% plot_od(od, filename, cmap_od, figure_width, figure_height)
%
% Plots the given OD map, saves the image if a filename is given, using
% the given colormap (default is gray), with figure widths and heights
% given in centimeters (defaults are 8 and 6).
%
% Requires export_fig:
% http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

if nargin < 2
    filename = nan;
end
if nargin < 3
    cmap_od = gray(17);
    cmap_od = cmap_od(1:16,:);
end
if nargin < 5
    figure_width = 8;
    figure_height = (3/4)*figure_width;
end

figure;
imagesc(od);
colormap(cmap_od);
axis equal;
axis tight;
axis off;
set(gca, 'YDir', 'Reverse');
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0.1, 0.1, -0.1, -0.1]);
set(gcf, 'Color', 'w', 'Units', 'Centimeters', 'Position', [10 10 figure_width figure_height]);
if ~isnan(filename)
    export_fig(filename);
    close;
end
