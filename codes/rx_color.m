function [color] = rx_color(rr)
% return legend color for plots for any receiver
%INPUT: receiver GRID ID
%OUTPUT: color
    switch rr 
        case 'grid108'  
            color = [0 0 1]; %blue
        case 'grid112'
            color = [0 1 0]; %green
        case 'grid154'
            color = [0.5 0.5 0.5]; %gray
        case 'grid160'
            color = [0 0 0]; %black
        case 'grid161'
            color = [0.1 0.6 0.3]; %dark green
        case 'grid162'
            color = [0 1 1]; %cyan
        case 'grid163'
            color = [1 0 1]; %magenta
        case 'ASTRArx'
            color = [1 0 0]; %red
    end
end

