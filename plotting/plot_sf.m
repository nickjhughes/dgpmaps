
function plot_sf(sf, filename, cmap_sf, figure_width, figure_height)
%PLOT_SF Plot the given SF map.
%
% plot_sf(od, filename, cmap_od, figure_width, figure_height)
%
% Plots the given SF map, saves the image if a filename is given, using
% the given colormap (default is jet), with figure widths and heights
% given in centimeters (defaults are 8 and 6).
%
% Requires export_fig:
% http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

if nargin < 2
    filename = nan;
end
if nargin < 3
    cmap_sf = jet(16);
end
if nargin < 5
    figure_width = 8;
    figure_height = (3/4)*figure_width;
end

sfmap = sf;
sfmap = 1.5*((sfmap-min(sfmap(:)))/max(sfmap(:)-min(sfmap(:))))+0.5;
d = log10(sfmap);
mn = min(d(:));
rnge = max(d(:))-mn;
d = 1+63*(d-mn)/rnge;
imagesc(d, [1+63*(log10(0.5)-mn)/rnge, 1+63*(log10(2)-mn)/rnge]);
colormap(cmap_sf);
axis equal;
axis tight;
axis off;
set(gca, 'YDir', 'Reverse');
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0.1, 0.1, -0.1, -0.1]);
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Centimeters', 'Position', [10 10 figure_width figure_height]);
if ~isnan(filename)
    export_fig(filename);
    close;
end
