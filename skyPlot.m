function hpol = skyPlot(varargin) 

load skyplotdata_PRN27_2015_076_zoom1.mat
prn1 = prnprn(1);
az1 = azprn(1,:);
el1 = elprn(1,:);
load skyplotdata_PRN22_2015_076_zoom1.mat
prn2 = prnprn(1,:);
az2 = azprn(1,:);
el2 = elprn(1,:);
load skyplotdata_PRN18_2015_076_zoom1.mat
prn3 = prnprn(1,:);
az3 = azprn(1,:);
el3 = elprn(1,:);
prn = [beamid;prn1;prn2;prn3];
% az = [beamAZ;[NaN az1];az2;az3];
az = [[NaN(4,1) beamAZ];az1;az2;[NaN az3]];
% el = [beamEL;[NaN el1];el2;el3];
el = [[NaN(4,1) beamEL];el1;el2;[NaN el3]];
varargin = {az,el,prn};

% load skyplotdata_PRN29_2014_051.mat
% prn = prnprn(1);
% az = azprn(1,:);
% el = elprn(1,:);
% prn = [beamid;prn];
% az = [beamAZ;az];
% el = [beamEL;el];
% 
% varargin = {az,el,prn};

% AZp = [repmat(AZb,[1 length(telist)]);AZ'*180/pi];
% ELp = [repmat(ELb,[1 length(telist)]);EL'*180/pi];
% prnp = [beamid;prn*ones(size(rcvr_op,1),1)];
% varargin = {az,el,prn};

%Function plots "sky view" from the receiver perspective.  
% 
% h = skyPlot(AZ, EL, PRN, line_style) 
% 
%   Inputs: 
%       AZ              - contains satellite azimuth angles. It is a 2D 
%                       matrix. One line contains data of one satellite. 
%                       The columns are the calculated azimuth values. 
%       EL              - contains satellite elevation angles. It is a 2D 
%                       matrix. One line contains data of one satellite. 
%                       The columns are the calculated elevations. 
%       PRN             - a row vector containing PRN numbers of the 
%                       satellites. 
%       line_style      - line style of the plot. The same style will be 
%                       used to plot all satellite positions (including 
%                       color).  
%   Outputs: 
%       h               - handle to the plot 
% 7 Aug 2013 S. Datta-Barua modified original at L171 to print prn text at 
% locations at the last point of the corresponding prn.

%-------------------------------------------------------------------------- 
%                           SoftGNSS v3.0 
%  
% Copyright (C) Darius Plausinaitis and Kristin Larson 
% Written by Darius Plausinaitis and Kristin Larson 
%-------------------------------------------------------------------------- 
%This program is free software; you can redistribute it and/or 
%modify it under the terms of the GNU General Public License 
%as published by the Free Software Foundation; either version 2 
%of the License, or (at your option) any later version. 
% 
%This program is distributed in the hope that it will be useful, 
%but WITHOUT ANY WARRANTY; without even the implied warranty of 
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
%GNU General Public License for more details. 
% 
%You should have received a copy of the GNU General Public License 
%along with this program; if not, write to the Free Software 
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, 
%USA. 
%-------------------------------------------------------------------------- 
 
%CVS record: 
%$Id: skyPlot.m,v 1.1.2.5 2006/08/18 11:41:57 dpl Exp $ 
 
%% Check arguments and sort them ========================================== 
[hAxis, args, nargs] = axescheck(varargin{:}); 
 
if nargs < 3 || nargs > 4 
    error('Requires 3 or 4 data arguments.') 
elseif nargs == 3 
    [az, el, prn]   = deal(args{1:3}); 
    line_style      = 'auto';     
else 
    [az, el, prn, line_style] = deal(args{1:4}); 
end 
 
if ischar(az) || ischar(el) || ischar(prn) 
    error('AZ and EL must be numeric.'); 
end 
 
if ~isequal(size(az), size(el)) 
    error('AZ and EL must be same size.'); 
end 
 
%% Prepare axis =========================================================== 
hAxis = newplot(hAxis); 
 
%--- Get x-axis text color so grid is in same color ----------------------- 
tc = get(hAxis, 'xcolor'); 
 
hold(hAxis, 'on'); 
 
%--- Plot white background ------------------------------------------------ 
rectangle('position', [-90, -90, 180, 180], ... 
          'Curvature', [1 1], ... 
          'facecolor', 'white', ... 
          'edgecolor', tc); 
 
%% Plot spokes ============================================================ 
 
%--- Find spoke angles ---------------------------------------------------- 
% Only 6 lines are needed to divide circle into 12 parts 
th = (1:6) * 2*pi / 12; 
 
%--- Convert spoke end point coordinate to Cartesian system --------------- 
cst = cos(th); snt = sin(th); 
cs = [cst; -cst]; 
sn = [snt; -snt]; 
 
%--- Plot the spoke lines ------------------------------------------------- 
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ... 
    'handlevisibility', 'off'); 
 
%% Annotate spokes in degrees ============================================= 
rt = 1.1 * 90; 
 
for i = 1:max(size(th)) 
 
    %--- Write text in the first half of the plot ------------------------- 
    text(rt*snt(i), rt*cst(i), int2str(i*30), ... 
        'horizontalalignment', 'center', 'handlevisibility', 'off', 'FontSize',14); 
 
    if i == max(size(th)) 
        loc = int2str(0); 
    else 
        loc = int2str(180 + i*30); 
    end 
 
    %--- Write text in the opposite half of the plot ---------------------- 
    text(-rt*snt(i), -rt*cst(i), loc, ... 
        'handlevisibility', 'off', 'horizontalalignment', 'center', 'FontSize',14); 
end 
 
%% Plot elevation grid ==================================================== 
 
%--- Define a "unit" radius circle ---------------------------------------- 
th = 0 : pi/50 : 2*pi; 
xunit = cos(th); 
yunit = sin(th); 
 
%--- Plot elevation grid lines and tick text ------------------------------ 
for elevation = 0 : 15 : 90 
    elevationSpherical = 90*cos((pi/180) * elevation); 
 
    line(yunit * elevationSpherical, xunit * elevationSpherical, ... 
        'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ... 
        'handlevisibility', 'off'); 
 
    text(0, elevationSpherical, num2str(elevation), ... 
        'BackgroundColor', 'white', 'horizontalalignment','center', ... 
        'handlevisibility', 'off', 'Fontsize',14); 
end 
 
%--- Set view to 2-D ------------------------------------------------------ 
view(0, 90); 
 
%--- Set axis limits ------------------------------------------------------ 
%save some space for the title 
axis([-95 95 -90 101]); 
 
%% Transform elevation angle to a distance to the center of the plot ------ 
elSpherical = 90*cos(el * pi/180); 
 
%--- Transform data to Cartesian coordinates ------------------------------ 
yy = elSpherical .* cos(az * pi/180); 
xx = elSpherical .* sin(az * pi/180); 
 
%% Plot data on top of the grid =========================================== 
for i = 1:length(prn)
    switch prn(i)
        case {23,27}
            color = [1, 0, 1];
            marker = 'd';
        case {13,22}
            color = [1, 0, 0];
            marker = '^';
        case {10,18}
            color = [0, 0, 0];
            marker = 'o';
        otherwise
            color =  [0, 0.447, 0.741];
            marker = 's';
    end
    if strcmp(line_style, 'auto') 
        %--- Plot with "default" line style ----------------------------------- 
        hpol = plot(hAxis, xx(i,:)', yy(i,:)', 'color', color); 
    else 
        %--- Plot with user specified line style ------------------------------ 
        % The same line style and color will be used for all satellites 
        hpol = plot(hAxis, xx(i,:)', yy(i,:)', line_style); 
    end 
 
    %--- Mark the initial position of the satellite ------------------------------ 
%     plot(hAxis, xx(i,1)', yy(i,1)',marker,'MarkerSize', 9, 'MarkerFaceColor', color); 

    %--- Mark the last position of the satellite ------------------------------ 
    plot(hAxis, xx(i,end)', yy(i,end)', ...
        marker ,'MarkerSize', 9, ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'none');
end

%--- Place satellite PRN numbers at the latest position ------------------- 
map1=colormap;map2=resample(map1,length(prn),size(map1,1));map2(find(map2<0))=0;map2(find(map2>1))=1;
for i = 1:length(prn) 
    if(prn(i) ~= 0 && prn(i) <= 32) 
        % The empthy space is used to place the text a side of the last 
        % point. This solution results in constant offset even if a zoom 
        % is used. 
%         text(xx(end, end), yy(end, end), ['  ', int2str(prn(i))], 'color',map2(i,:)); 
%         text(xx(i, end), yy(i, end), ['  ', int2str(prn(i))], 'color',map2(i,:),'FontSize',14, 'FontWeight','bold'); 
        text(xx(i, end), yy(i, end), ['  ', int2str(prn(i))]); 
    end 
end 
 
%--- Make sure both axis have the same data aspect ratio ------------------ 
axis(hAxis, 'equal'); 
 
%--- Switch off the standard Cartesian axis ------------------------------- 
axis(hAxis, 'off'); 

[~, op_path] = ver_chk;
% saveas(gcf,[op_path, '/thesis/figs/', ...
%     'PRN27_22_18_2015_076_zoom1_0-1800s_after_1100UT_skyplot.eps'],'epsc2');