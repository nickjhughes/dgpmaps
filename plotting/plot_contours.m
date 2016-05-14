
function plot_contours(op, od, filename, nangles, cmap_op_contours, col_od, pinwheels, figure_width, figure_height)
% Plot the contours of the given OP and OD maps. If filename is given, save to file.

if nargin < 3
    filename = nan;
end
if nargin < 4
    nangles = 8;
end
if nargin < 5
    cmap_op_contours = hsv(nangles);
end
if nargin < 6
    col_od = 'k';
end
if nargin < 7
    pinwheels = locate_pinwheels(op);
end
if nargin < 9
    figure_width = 8;
    figure_height = (3/4)*figure_width;
end

angles = linspace(-90, 90, nangles+1);
angles = angles(1:end-1);
c = op_contours(op, angles);
figure;
hold on;
cmap = cmap_op_contours;
for j = 1:length(c)
    plot(c{j}(1,:), c{j}(2,:), 'Color', cmap(j,:), 'LineWidth', 1);
end
c = contourc(od, [0 0]);
c(:,c(1,:)==0) = nan;
plot(c(1,:), c(2,:), 'Color', col_od, 'LineWidth', 2);
if ~isempty(pinwheels)
    plot(pinwheels(:,1), pinwheels(:,2), '.', 'Color', col_od, 'MarkerSize', 20);
end
hold off;
set(gca, 'YDir', 'Reverse');
axis equal;
axis tight;
axis off;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'Centimeters', 'Position', [10, 10, figure_width, figure_height]);
if ~isnan(filename)
    export_fig(filename);
    close;
end
