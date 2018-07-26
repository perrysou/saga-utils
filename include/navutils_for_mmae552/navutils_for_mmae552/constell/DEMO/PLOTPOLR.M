function hpol = plotpolr(az,el,line_style,prn_label,label_color,spoke_label)

% hpol = plotpolr(az,el,line_style,prn_label,label_color,spoke_label);
%
% Azimuth-Elevation polar coordinate plot. Plots a line based on azimuth and 
% elevation data with optional inputs for labeling and coloring of the line.
%
% Inputs:
%   az           - azimuth vector (rad) (nx1)
%                   valid values are -2*pi -> 2*pi
%   el           - elevation vector (rad) (nx1)
%                   valid values are -pi/2 -> pi/2
%   line_style   - line style to use for this line (string) (optional). 
%                   These are the Matlab character designations for line styles 
%                   [. o x + - * : -. --], default is the Matlab 'auto'
%   prn_label    - label to identify the az/el line (string) (optional)
%                   default is no label
%   label_color  - color for the line and label (string) (optional).  
%                   These are the Matlab letter designations for colors 
%                   [y m c r g b w k], default is the Matlab 'auto'
%   spoke_label  - labels for the North, East, South,and West spokes of the 
%                   plot. (4xn string) (optional).  
%                   default = str2mat('N','E','S','W')
%                  
% Outputs:
%   hpol         - handle to the plot created
%
% See also PLOTAZEL, NED2AZEL, ECEF2NED

% Written by: Jimmy LaMance 11/7/96 
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
p_handle=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,6,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PLOTPOLR for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check if the line style is provided, if not use the default of auto
if nargin < 3            % no line style variable is provided
  line_style = 'auto';   % set to auto mode (default)
end % if nargin < 3
  
% check if the prn label is provided, if not use the default of '' (no label)
if nargin < 4            % no prn label variable is provided
  prn_label = '';        % set to '' (default)
end % if nargin < 4

% check if the label color is provided, if not use the default of auto
if nargin < 5            % no line style variable is provided
  label_color = '';      % set to auto mode (default)
end % if nargin < 5 

% check that the input spoke labels are correct dimension or set to the 
% default values
if nargin < 6
  spoke_label = str2mat('N','E','S','W');
end % if nargin >= 6          

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'PLOTPOLR';

  % Develop the error checking structure with required dimension, matching
  % dinemsion flags, and input dinemsions.
  estruct.variable(1).name = 'az';
  estruct.variable(1).req_dim = [901 1];
  estruct.variable(1).var = az;
  estruct.variable(1).type = 'ANGLE_RAD';
  
  estruct.variable(2).name = 'el';
  estruct.variable(2).req_dim = [901 1];
  estruct.variable(2).var = el;
  estruct.variable(2).type = 'ELEVATION_RAD';

  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% check the current version of Matlab in use
matlab_version = version;      % full version string

% get the first digit (either 4 or 5), or report an error
version_number = str2num(matlab_version(1));

% check that this is version 4 or 5, if not report not supported and return
if version_number ~= 4 & version_number ~= 5
  fprintf('This function is only supported under Matlab versions 4 and 5.\n') ;
  return
end % if version_number ~= 4 & version_number ~= 5

% get current axis and the hold state
cax = gca;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% rotate x (az) 90 degrees to line up with the label defined above
az = -az + pi/2;

% convert elevation to degrees (azimuth will stay in radians)
el = el * 180 / pi;

% convert elevation to declination
el = 90 - el;    

% transform data to Cartesian coordinates.
xx = el.*cos(az);
yy = el.*sin(az);

% black out all the data that is not on the grid
% this section of code is here so that the grid is drawn 
% after the points are plotted and blacked out
I_out_plot = find(sqrt(xx.^2 + yy.^2) >= 90);

% if there are points off the plot, black them out (radius of 90 degree in el)
% this is done by setting the polt location to inf which does not show
% on a plot, a nice feature in Matlab
if any(I_out_plot)    
  xx(I_out_plot) = ones(size(I_out_plot,1),1) * inf;
  xy(I_out_plot) = ones(size(I_out_plot,1),1) * inf;;
end % if any(I_out_plot)

% find the last data point
last = size(xx,1);

% put the x/y data pairs on the plot based on user supplied line styles
if strcmp(line_style,'auto')
    
  hpol = plot(xx,yy);              % plot command with auto line style
    
  if strcmp(label_color,'')     % no color provided, let Matlab choose
    % put the end point as an x with a Matlab chosen color
    qte = text(xx(last),yy(last),'x');  
  else                          % color provided, user the user supplied color
    % put the end point as an x with the desired color
    qte = text(xx(last),yy(last),'x','Color',label_color);   
  end % if strcmp(label_color,'')

else     % line style is provided

  if strcmp(label_color,'')     % no color provided, let Matlab choose
    % generate the plot and put the end point as an x with a Matlab chosen color
    if length(xx) > 1
      hpol = plot(xx,yy,line_style);
      qte = text(xx(last),yy(last),'x');    
    end % if length(xx) > 1
  else                          % color provided, user the user supplied color
    % generate the plot and put the end point as an x with the desired color
    if length(xx) > 1
      hpol = plot(xx,yy,line_style,'Color',label_color);
      qte = text(xx(last),yy(last),'x','Color',label_color);   
    end % if length(xx) > 1
 end % if strcmp(label_color,'')

end % if strcmp(line_style,'auto')

% set the background color to black
set(gcf,'Color','k');

% update the current axies
cax = gca;

% set x-axis and y-axis text color to white 
set(cax,'Color','k','XColor','w','YColor','w');

% set the variable tc (used for spoke color) to white to match the axis colors
tc = 'w';

% end point symbols were placed on the last point in the data
% if the end point is determined to be off the plot, then...
if any(I_out_plot)
  I_last = find(I_out_plot == last);
  if any(I_last)                  % need to clear an end symbol
    qte = text(xx(last),yy(last),'x','Color','k');  % end point symbol is an x
  end % if any(I_last)
end % if any(I_out_plot)

% label the line with a satellite number, 
% put the prn number in the middle of the arc
I_in_plot = find(sqrt(xx.^2 + yy.^2) < 90);

% verify that there is a place to put the label and that a 
% valid label is available
if any(I_in_plot) & ~strcmp(prn_label,'')
  index = fix(size(I_in_plot,1) / 2);    % index to point to add label to
  
  if index == 0     % just in case there is only one point that passed
    index = 1; 
  end % if index == 0
               
  % compute location for the label
  x_label_point = xx(I_in_plot(index));
  y_label_point = yy(I_in_plot(index));
  
  qts = text(x_label_point,y_label_point,prn_label,'Color',label_color);         % put the label on
  
  % clear out unneeded variables
  clear I_in_plot index         
  
end % if any(I_in_plot)  

% update the current axies
cax = gca;

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
  hold on;              
  max_angle = 2*pi;      % 0-360 degrees
  max_radius = 90;       % 0-90 degrees
  hhh=plot([0 max_angle],[0 max_radius],'w');
  v = [get(cax,'xlim') get(cax,'ylim')];
  ticks = length(get(cax,'ytick'));
  delete(hhh);
  
  % check radial limits and ticks
  rmin = 0;    
  rmax = max_radius; 
  rticks = 3;   % force only three rings

% define a circle with pi/50 azmuthal resolution
  th = 0:pi/50:2*pi;
  xunit = cos(th);
  yunit = sin(th);

% now really force points on x/y axes to lie on them exactly
  inds = [1:(length(th)-1)/4:length(th)];
  xunits(inds(2:2:4)) = zeros(2,1);
  yunits(inds(1:2:5)) = zeros(3,1);

  rinc = (rmax-rmin)/rticks;
  for i=(rmin+rinc):rinc:rmax      % control of elevation
    plot(xunit*i,yunit*i,'-','color',tc,'linewidth',1);
    text(0,i+rinc/20,['  ' num2str(abs(i - 90))],  ...
         'verticalalignment','bottom','Color','w' );
  end % for i=(rmin+rinc):rinc:rmax

% plot spokes
  th = (1:6)*2*pi/12;
  th = (1:6)*2*pi/12;
  cst = sin(th); snt = cos(th);    % formualtion of, 0 - N, 90 - E, etc
  cs = [-cst; cst];
  sn = [-snt; snt];
  plot(rmax*cs,rmax*sn,'-','color',tc,'linewidth',1);

% annotate spokes in degrees (except the N S E and W directions)
  rt = 1.1*rmax;
  
  for i = 1:max(size(th))
    if rem(th(i), pi/2) ~= 0
      text(rt*cst(i),rt*snt(i),int2str(i*30), ...
           'horizontalalignment','center', ...
           'EraseMode','none','Color','w');
      if i == max(size(th))
        loc = int2str(0);
      else
        loc = int2str(180+i*30);
      end % if i == max(size(th))
      
      text(-rt*cst(i),-rt*snt(i),loc, ...
           'horizontalalignment','center', ...
           'EraseMode','none','Color','w');
    end % if rem(th(i), pi/2) ~= 0
  end % for i = 1:max(size(th))

% annotate the N,S,E, and W spokes using the string supplied in spoke_label
  th = (0:3)*2*pi/4;
  cst = sin(th); snt = cos(th);    
  
  rt = 1.15*rmax;

  for i = 1:max(size(th))
    text(rt*cst(i),rt*snt(i),spoke_label(i,:), ...
         'horizontalalignment','center', ...
         'EraseMode','none','Color','w');
  end % for i = 1:max(size(th))

  % add a notation that states that the end points are labeled with an x
  gt1 = text(-1.75 * rmax, -1.15 * rmax, '''x'' denotes last point in time');
  gt2 = text(-1.85 * rmax, -1.25 * rmax, 'elevation is 90 degrees at center'); 
  set(gt1,'Color','w');
  set(gt2,'Color','w');

  % set the axis color to white
  set(cax,'color','k');
  set(cax,'XColor','w','YColor','w');
  
  % set viewto 2-D
  view(0,90);

% code from previous version, not needed in the future???
%  if version_number < 5
%    axis(rmax*[-1 1 -1.1 1.1]);
%    axis(rmax*[-1 1 -.9 .9]);
%  end % if version_number < 5

end % if ~hold_state

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight );

% reset the axis if not in a hold state
if ~hold_state
  % set the x and y axis to equal so that a circle is formed
  axis('equal');
  
  % force the axis limits (required to show correctly in Matlab 5)
  set(cax,'XLim',[-90 90],'YLim',[-105 110]);          
  
  % hide the axes
  axis('off');
end % if ~hold_state

% reset hold state
if ~hold_state 
  set(cax,'NextPlot',next); 
end % if ~hold_state

% return the handle to the current figure
hpol = gcf;

%%%%% END ALGORITHM CODE %%%%%

% end PLOTPOLR
