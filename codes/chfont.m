function chfont(scale,type)
% a oneliner that changes the font size of the current fig
switch nargin    
    case 1
        set(findall(gcf,'-property','FontSize'),'FontSize',scale);
    case 2
        set(findall(gcf,'-property','FontSize'),'FontSize',scale);
        set(findall(gcf,'-property','FontSize'),'FontName',type);
end
end

