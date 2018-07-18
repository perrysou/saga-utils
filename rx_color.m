function [color] = rx_color(rr)
% return legend color for plots for any receiver
%INPUT: receiver GRID ID
%OUTPUT: color
switch rr
    case {'grid108', 'IIT-1'}
        color = [0, 0, 1]; %blue
    case {'grid112'}
        color = [0, 1, 0]; %green
    case {'grid154', 'IIT-16'}
        color = [0.5, 0.5, 0.5]; %gray
    case {'grid160', 'IIT-9'}
        color = [0, 0, 0]; %black
    case {'grid161', 'IIT-15'}
        color = [0.1, 0.6, 0.3]; %dark green
    case {'grid162', 'IIT-11'}
        color = [0, 1, 1]; %cyan
    case {'grid163', 'IIT-3'}
        color = [1, 0, 1]; %magenta
    case {'ASTRArx', 'IIT-13'}
        color = [1, 0, 0]; %red
end
end
