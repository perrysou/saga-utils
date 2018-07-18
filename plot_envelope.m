function [ax] = plot_envelope(ax, x, y, e, c)

x = x(~isnan(y));
e = e(~isnan(y));
y = y(~isnan(y));
% errorbar(ax,x,y,e,'.');

hold on;
if isrow(x) && isrow(y) && isrow(e)
    x = x';
    y = y';
    e = e';
end
fill([x; flipud(x)], [y - e; flipud(y+e)], ...
    c, 'facecolor', c, 'facealpha', 0.1, 'edgecolor', 'none');
% patch([x; flip(x)], [y; flip(y+e)], c, 'facecolor', c, 'facealpha', 0.1, 'edgecolor', 'none');
