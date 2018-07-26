function current_fig = orb_anim(alm,start_time,animation_dt,sta_loc,...
                                 sta_name,mask)

% current_fig = orb_anim(alm,start_time,animation_dt,sta_loc,sta_name,mask);
%
% Animate the orbits and ground stations over a Mercator map of the Earth.
%
% Inputs:
%   alm          - almanac (nx13) or ephemeris (nx22) data for 
%                   satellites to be animated
%   start_time   - GPS start time, [GPS_weeks GPS_seconds] (1x2)
%                   valid GPS_week values are 1-3640 (years 1980-2050)
%                   valid GPS_sec values are 0-604799
%   animation_dt - time step for the animation (1x1) (sec) (optional)
%                   default = 60 
%   sta_loc      - Earth station locations [lat, lon, hgt] (nx3) (optional)
%                   lat and lon (East) are in rad, hgt is in meters
%   sta_name     - Earth station names (nxm) (string) (optional)
%   mask         - station visiblity mask (rad) (optional) default = 0.0 
%                   mask can be 1x1, nx4, nx5 in the form 
%                   [sta_num min_el min_az max_az]   (nx4)
%                   [sta_num min_el max_el min_az max_az]   (nx5)
%                   with n elevation/azimuth triples, 
%                   where sta_num corresponds to the sta_loc/sta_name index.
%                   if overlapping visibility masks are given, the least 
%                   stringent mask will be applied
% Output:
%   current_fig  - handle to the resulting figure
%
% See also PLOTSKY, WRITEMPG 

% Written by: Jimmy LaMance 4/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, ALM2GEPH, GPS2UTC, PROPGEPH, ECEF2LLA,
%                   LLA2ECEF, MASK_STA, GUISLIDR, LOS, ECEF2NED, 
%                   NED2AZEL, VIS_DATA

% WGS-84 constants
RADIUS_EARTH = 6378137.0;        % radius of the Earth WGS-84 value

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
current_fig=[];

% set the animation time step
if nargin < 3
  animation_time_step = 60; % set animation time step to the default of 60 minutes
  animation_dt = 60;
else
  animation_time_step = animation_dt;   % use the input value  
end % if nargin < 3

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,6,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ORB_ANIM for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check the size of the station locations, if provided
if nargin < 5      % station location provided
  % set the flag for the use of base stations to false
  base_station_loc_flag = 0;
else 
  base_station_loc_flag = 1; 
  estruct.variable(5).name = 'sta_loc';
  estruct.variable(5).req_dim = [901 3];
  estruct.variable(5).var = sta_loc;
end % if nargin < 5

% check that the station names are strings, if provided
if nargin >= 5      % station names provided
  % set the flag for the use of base stations to true
  base_station_name_flag = 1;
else                 % no station names provided
  % set the flag for the use of base names to false
  base_station_name_flag = 0;
end % if nargin >= 5

% set the default station mask if not provided
% check that the size of the mask input is valid, if provided
if nargin < 6
  mask = 0;    % set mask to the default value of 0
end % if nargin < 6

estruct.func_name = 'ORB_ANIM';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'alm';
estruct.variable(1).req_dim = [903 13];
estruct.variable(1).var = alm;
  
estruct.variable(2).name = 'start_time';
estruct.variable(2).req_dim = [1 2];
estruct.variable(2).var = start_time;
estruct.variable(2).type = 'GPS_TIME';
  
estruct.variable(3).name = 'animation_dt';
estruct.variable(3).req_dim = [1 1];
estruct.variable(3).var = animation_dt;
  
estruct.variable(4).name = 'mask';
estruct.variable(4).req_dim = [1 1; 905 4; 905 5];
estruct.variable(4).var = mask;
  
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% check that the base station location longitudes are all positive
if base_station_loc_flag == 1 
  I_neg_lon = find(sta_loc(:,2) < 0);
  if any(I_neg_lon)
    sta_loc(I_neg_lon,2) = sta_loc(I_neg_lon,2) + 2 * pi;
  end % if any(I_neg_lon)
end % if base_station_loc_flag == 1

% set come default colors, this is done here for ease of modification
% colors are percentage of RGB
base_station_label_color = [0 0 1];
base_station_marker_color = [0 0 1];
base_station_mask_color = [.8 .8 .8];
visible_sat_color = [0 1 0];
non_visible_sat_color = [1 0 0];
topo_line_color = [0 0 0];

erase_type = 'xor';
                                                             
% convert the almanac to an ephemeris, if necessary
if size(alm,2) == 13
  ephem = alm2geph(alm);  
elseif size(alm,2) == 22
  ephem = alm;
else
  fprintf('Unknown almanac input format.\n');
  fprintf('This should have been detected in the code previously.\n');
  fprintf('Contact customer support.\n');
  return
end % if size(alm,2) == 13

% set the starting time
gps_time = start_time;
gps_time0 = gps_time;

% find the number of satellites to be animated
num_animated_satellites = size(alm,1);   % initialize the number of satellites

% create the string labels to be used for each vehicle based on the input prn 
% use 2, 3, or 4-digit satellite numbers depending on the maximum value 
% in the almanac
max_satellite_number = max(alm(:,1));

% Round the first column of the almanac to get integer satellite numbers
a_round = round(alm(:,1));
prn_diff = a_round - alm(:,1);
I_diff_non_zero = find(prn_diff ~= 0);
if ~isempty(I_diff_non_zero)
  alm(I_diff_non_zero,1) = a_round(I_diff_non_zero,1);
  fprintf('Non-integer satellite numbers have been rounded in ORB_ANIM.\n')
end % if ~isempty(I_diff_non_zero)

for i = 1:num_animated_satellites
  if max_satellite_number < 100
    head_labels(i,1:2) = sprintf('%2d',alm(i,1));
  elseif max_satellite_number < 1000
    head_labels(i,1:3) = sprintf('%3d',alm(i,1));
  else
    head_labels(i,1:4) = sprintf('%4d',alm(i,1));
  end % if max_satellite_number < 100
end

% set the output figure handle to the current figure
current_fig = gcf;

% generate the first time labels to put on the plot and set up 
% the properties of the text handles
% use erase mode xor for speed and smaller font size for ease of display
utc_time = gps2utc(gps_time);     % compute UTC from GPS time

% turn on the zoom property
zoom on                 

% set the background and x an y-axis colors 
set(gca,'Color','w','XColor','k','YColor','k','box','on');
  
% turn the grid off
set(gca,'XGrid','off','YGrid','off');
 
set(gca,'XLim',[0 360],'YLim',[-90 90]);

set(gcf,'Resize','off');
  
% set the figure color to white
set(gcf,'Color','w')
  
% get the handle to the title
title_handle = get(gca,'Title');
  
% set the title to be black
set(title_handle,'Color','k');

% set the renderer mode 
%set(gcf,'RendererMode','manual','Renderer','zbuffer');

% set the current dither map and set the dither map mode to manual 
% so that it does't go through the dither process each time a new 
% object is added to the figure
%d_map = get(gcf,'Dithermap');
%set(gcf,'DithermapMode','manual')

% set the interruptible property for the different versions
set(gcf,'Interruptible','Off'); 

% set animation properties
set(gcf,'BackingStore','On');  % set the figure backstoring to on

% set the pointer to be an arrow so that the user feels the window is active
set(gcf,'Pointer','Arrow');

% set draw mode to fast (speeds the animation)
set(gca,'DrawMode','fast'); 

%label the x & y-axis and the plot title
ylabel('Latitude (degrees)')
xlabel('Longitude (degrees)') 
title_text = sprintf('Satellite Animation %02d/%02d/%02d  %d:%02d:%02.0f ',...
                     utc_time(2),utc_time(3),utc_time(1)-1900,utc_time(4:6));
th = title(title_text); 
set(th,'EraseMode','xor');

% get the handle to the title text
th = get(gca,'title'); 

% set the title erase mode to the erase type being used
% set the interpreter (version 5) to none to speed animation
set(th,'Interpreter','none');

% set the position of the animation window and the GUI control window
% first get the original screen size and set the size to pixels
orig_units = get(0,'Units');
set(0,'units','pixels');
set(gcf,'units','pixels');     % set the figure units too

% get the screen size
screen_size = get(0,'ScreenSize'); 

% reset the base property value
set(0,'units',orig_units);

% set the width to be 3/4 of the screen
width = .75 * screen_size(3);

% use an aspect ration of 1.3 to establish the height
height = fix(width / 1.3);  

% set the position of the animation window to be in the top left corner
an_win_pos(1) = 10;                             % x start

an_win_pos(2) = screen_size(4) - height - 60;   % y start 

an_win_pos(3) = width;
an_win_pos(4) = height;

% set the tick mark locations and labels
set(gca,'Xtick',[0 60 120 180 240 300 360]);
set(gca,'Ytick',[-90 -60 -30 0 30 60 90]);

set(gca,'XtickLabel',str2mat('0','60','120','180','240','300','360'));
set(gca,'YtickLabel',str2mat('-90','-60','-30','0','30','60','90'));

% set position property to these values
set(gcf,'position',an_win_pos);

% name the figure 
set(gcf,'Name','Constellation Toolbox Animation Window');    
   
% turn off the automatic figure numbering
set(gcf,'numbertitle','off');

% compute the satellite positions for this time
[t_out, prn, x_ecef, v] = propgeph(ephem, gps_time);  
  
% covert ECEF positions to lla in degrees
r2d = 180 / pi;     % conversion from degrees to radians
x = ecef2lla(x_ecef) * r2d;

% do some computations on the base station location data for plotting
if base_station_loc_flag == 1
  % determine the number of tracking stations
  number_of_stations = size(sta_loc,1);

  % rename the station location with radians, keep the degrees version
  sta_loc_rad = sta_loc;                    

  % convert station location in degrees (that's the way the plot is done)
  sta_loc = [sta_loc_rad(:,1:2) * r2d sta_loc_rad(:,3)];
  
  % compute ECEF station locations
  sta_loc_ecef = lla2ecef(sta_loc_rad);

  % make sure the longitudes are all positive
  I = find(x(:,2) < 0);
  if any(I)
    x(I,2) = x(I,2) + 360;
    clear I
  end % if any(I)
end % if base_station_loc_flag == 1
    
hold on            

% initialize the drawing      
for n = 1:num_animated_satellites       
  head(n) = text(inf,inf,head_labels(n,:),'color',non_visible_sat_color,...
                  'erase',erase_type);  
%  head = text(inf,inf,head_labels,'color',non_visible_sat_color,...
%                  'erase',erase_type);  
  
end % for n         

% set the horizontal alignment of the satellite indicators to be centered
set(head,'HorizontalAlignment','center')

% draw the ground station mask on the Earth, starting with the elevation
% mean altitude of the satellites of interest, ephemeris(:,5) is the sqrt
% of the semi-major axis
orb_height_earth = mean(ephem(:,5).^2) - RADIUS_EARTH;     
    
% add the tracking stations if input
if base_station_loc_flag == 1
  if nargin > 3 
    for i = 1:number_of_stations;

      % find the masking elements that correspond to this station
      if length(mask) ~= 1              % if the nx4 version of masking is used
        I_mask = find(mask(:,1) == i);  % find the masking for this station
        
        if any(I_mask)     % if there is specific masking information
          this_mask = mask(I_mask,2:end);      % for this station use it
        else                                 % else
          this_mask = 0;                     % use the default masking os zero
        end % if any(I_mask)
          
      else                                   % else
        this_mask = [mask 0 2*pi];           % use the elevation mask only
      end % if size(mask,1) ~= 1

      % compute the elevation mask to draw on the ground
      [lat, lon] = mask_sta(sta_loc_rad(i,:),this_mask,orb_height_earth,0.001); 
    
      % make sure all of the lons are positive
      I = find(lon < 0);
      if any(I)
        lon(I) = lon(I) + 2 * pi;
      end % if any(I)
    
      % convert the lat and lon to degrees for plotting
      r2d = 180 / pi;
      lon = lon * r2d;
      lat = lat * r2d;                     

      % plot the ground station visiblity coverage masks
%      h_mask(i) = line(lon, lat,...
%                       'color',base_station_mask_color,...
%                       'markersize',5,'erasemode','none',...
%                       'LineStyle','-'); 
      h_mask = line(lon, lat,'marker','.',...
                       'color',base_station_marker_color,...
                       'markersize',5,'erasemode','none',...
                       'LineStyle','none'); 
%keyboard
%      h_mask(i) = patch(lon, lat,base_station_mask_color,...
%                        'facecolor',base_station_mask_color,...
%                        'edgecolor',base_station_label_color,...
%                        'erasemode','none'); 
      
      h_sta(i) = line(sta_loc(i,2),sta_loc(i,1));
      
      set(h_sta(i),'marker','+','color',base_station_marker_color,...
                  'markersize',5,'erasemode','none'); 
            
      if base_station_name_flag == 1  
        h_sta_text(i) = text(sta_loc(i,2),sta_loc(i,1),sta_name(i,:));
        set(h_sta_text(i), 'VerticalAlignment','bottom','fontsize',10,...
                           'erasemode','none','color',base_station_label_color); 
      end % if base_station_name_flag == 1

    end % for i = 1:size(sta_loc,1);
  end % if nargin > 3  
end % if base_station_loc_flag == 1

% load the world topo data
load topo                                

% generate the contour plot at sea level, plot the topo contour at 0
[c, topo_handles] = contour(1:360,-89:90,topo,[0 0],'k-');  
  
% set the erase property to none
set(topo_handles,'EraseMode','none','color',topo_line_color)

% resize the axes so we can put the animation control on the active figure
% window.  this is essential for proper opration in Matlab 5.  otherwise
% the animation is REALLY choppy because the inactive window does not get
% updated quickly enough   

set(gcf,'units','pixels');
fig_pos = get(gcf,'Position');
x_extent = fig_pos(3) - fig_pos(1);
y_extent = fig_pos(4) - fig_pos(2);
set(gca,'units','pixels');

fig_pos = [fix(x_extent*.1) fix(y_extent*.25) fix(x_extent*.75) fix(y_extent*.85)];
set(gca,'Position',fig_pos);

% set the user data values to zero
adjust_values = zeros(11,1);    

% open the GUI window for controlling the animation
hs = guislidr(gcf);     

% set the animation speed to the input value
adjust_values = get(hs,'userdata');
adjust_values(4) = animation_dt;            

% set the properties in the gui slider properties
set(hs,'userdata',adjust_values);

% set the initial speed value
guislidr(gcf,'SetSpeed');

% now set the position for the animation control window (bottom right)
% first get the size of the animation control window in pixels 
set (hs,'units','pixels')
ctrl_size_org = get(hs,'position');
ctrl_width = ctrl_size_org(3);
ctrl_height = ctrl_size_org(4); 

% compute the new window positions
ctrl_win_pos(1) = screen_size(3) - ctrl_width - 40;
ctrl_win_pos(2) = 35;
ctrl_win_pos(3) = ctrl_width;
ctrl_win_pos(4) = ctrl_height;

% put all of this in a while loop that is closed with the Done button
done_flag = 0;     
pause_flag = 0;
play_flag = 1; 

% bring the animation window to the front
figure(gcf)

drawnow

while done_flag == 0;  

  % compute the satellite positions for this time
  [t_sats, prn, x_ecef, v] = propgeph(ephem, gps_time);  
  
  % convert to lla
  x = ecef2lla(x_ecef);
  x(:,1:2) = x(:,1:2) * r2d;

  % make sure the longitudes are all positive
  I = find(x(:,2) < 0);
  if any(I)
    x(I,2) = x(I,2) + 360;
    clear I
  end % if any(I)
    
  % plot the head points by changing the position property
  pn = {'position'};
  pv = num2cell([x(:,2),x(:,1)],2);
  set(head,pn,pv);
  
%keyboard
  % put the time tag on the plot
  utc_time = gps2utc(gps_time);    % update the date and times

  % update the text in the title
  title_text = sprintf('Satellite Animation %d/%d/%d  %d:%02d:%02.0f ',...
                       utc_time(2),utc_time(3),utc_time(1)-1900,utc_time(4:6));
  set(th,'string',title_text);

  % compute LOS to each base station and change the color of the satellite
  % if visible from any station (if station data is input) 
  if base_station_loc_flag == 1
    I_pass = [];
    for ns = 1:number_of_stations 
      % compute LOS in ECEF coordinates
      [t_los,los_ecef,sta_num_temp,veh_num] = ...
                         los(t_sats(1,:),sta_loc_ecef(ns,:),t_sats,[prn x_ecef]);
    
      % convert to NED
      [los_ned] = ecef2ned(los_ecef, sta_loc_rad(ns,:));
    
      % compute az and el in radians
      [los_az, los_el] = ned2azel(los_ned);
    
      % find the masking elements that correspond to this station
      % if the nx4 version of masking is used
      if length(mask) ~= 1                 
        % find the masking for this station
        I_mask = find(mask(:,1) == ns);       
        
        % if there is specific masking information for this station
        if any(I_mask)                       
          this_mask = mask(I_mask,2:end);    % use it 
        else                                 % else
          this_mask = 0;                     % use the default masking of zero
        end % if any(I_mask)
          
      else                                   % else
        this_mask = mask;                    % use the elevation mask only
      end % if size(mask,1) ~= 1
      
      % find which satellites are visible
      [pass_az_el, I_pass_this_sta] = vis_data(this_mask,[los_az los_el]);
      I_pass = [I_pass; I_pass_this_sta];

    end

    % define a matrix with ones in the visible satellite location
    vis_matrix = zeros(num_animated_satellites,1);
    vis_matrix(I_pass) = ones(size(I_pass));
  
    clear I_pass
  
    I_pass = find(vis_matrix == 1);
    I_fail = find(vis_matrix == 0);
  
    % set all of the colors
    if ~isempty(I_fail)
      set(head(I_fail),'color',non_visible_sat_color);
    end % if any(I_fail)
  
    % change the color for the satellites in view
    if ~isempty(I_pass)
      set(head(I_pass),'color',visible_sat_color);
    end  

  end % if base_station_loc_flag == 1

  drawnow
  
  % check the value for the GUI control to adjust the speed
  % the userdata has the following meaning
  % column 4 - speed, 5 - play, 6 - pause, 7 done
  adjust_values = get(hs,'userdata'); 
  
  % check and adjust the speed
  if adjust_values(4) == 0
    adjust_values(4) = 1;
  end
  
  animation_dt = adjust_values(4);     % adjust(4) == speed
    
  % see if the pause flag has been set
  while adjust_values(6) ~= 0          % adjust(6) == pause
    pause(.1)
    adjust_values = get(hs,'userdata'); 
  end % while

  % see if the done flag has been set
  if adjust_values(7) ~= 0             % adjust(7) == done
    break
  end % if done_flag ~= 0
  
  % go to the next time step
  gps_time(2) = gps_time(2) + animation_dt;
  
  % check for week rollover
  if gps_time(2) >= 86400 * 7
    gps_time(1) = gps_time(1) + 1;
    gps_time(2) = gps_time(2) - 86400 * 7;
  end % if gps_time(2) >= 86400 * 7

  % see if the stop flag has been set
  if adjust_values(8) ~= 0             % adjust(8) == restart
    
    gps_time = gps_time0;   
    adjust_values(8) = 0;
    set(hs,'userdata',adjust_values);
    
  end % if done_flag ~= 0 
  
end % while done_flag == 0

hold off   % turn the figure hold off and we're done

close(current_fig);

%%%%% END ALGORITHM CODE %%%%%

% end ORB_ANIM     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  GUISLIDR function within ORB_ANIM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_fig = guislidr(h_fig,command_str)

% h_fig = guislidr(h_fig,command_str);
%
% This function displays and controls the animation control GUI.  
% Used as a support function for ORB_ANIM.
%
% Input:
%   h_fig       - handle to the figure window to locate the GUI
%   command_str - command string
% Output:
%   h_fig       - handle to the GUI figure window
%
% See also ORB_ANIM

% Written by: Jimmy LaMance 1/15/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% no error checking required in this function

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

if nargin < 2
    command_str = 'initialize';
end % if nargin < 2

if strcmp(command_str,'initialize')     % initialize the GUI window

  % check the current version of Matlab in use
  matlab_version = version;      % full version string

  % get the first digit (either 4 or 5), or report an error
  version_number = str2num(matlab_version(1));

  % check that this is version 4 or 5, if not report not supported and return
  if version_number ~= 4 & version_number ~= 5
    fprintf('This function is only supported under Matlab versions 4 and 5.\n') ;
    return
  end % if version_number ~= 4 & version_number ~= 5
        
  % set the units for placing the UIcontrols in normalized.  This will
  % make it a lot easier moving to and from different screen sizes
  
  % get the position of the current axis
  ax_pos = get(gca,'position');
  
  % check the current version of Matlab in use
  matlab_version = version;      % full version string

  % get the first digit (either 4 or 5), or report an error
  version_number = str2num(matlab_version(1));

  % check that this is version 4 or 5, if not report not supported and return
  if version_number ~= 4 & version_number ~= 5
    fprintf('This function is only supported under Matlab versions 4 and 5.\n') ;
    return
  end % if version_number ~= 4 & version_number ~= 5
  
  % set the width and length of the push buttons in pixels 
  % the buttons should be sized based on the Matlab version in use because
  % of the way the UI interface differes from version 4 to 5
  if version_number > 4
    button_width = 45;
    button_height = 20;
    bvs = 30;                % button_vertical_separation (easier to code)
  else
    button_width = 60;
    button_height = 25;
    bvs = 35;                % button_vertical_separation (easier to code)
  end % if version_number > 4
    
  % determine the offsets from the current axis to where the buttons will go
  x_offset = ax_pos(1) + ax_pos(3) + 20;
  y_offset = ax_pos(2) + ax_pos(4) - 60;
  
  % set the positions for the buttons
  h_play_pos = [x_offset y_offset button_width button_height];
  h_pause_pos = [x_offset y_offset-bvs button_width button_height];
  h_done_pos = [x_offset y_offset-2*bvs button_width button_height];
  h_stop_pos = [x_offset y_offset-3*bvs button_width button_height];
  
  % set the positions for the speed slider
  h_spdsld_pos = [35 5 170 15];   
  h_spdmin_pos = [10 20 55 20];
  h_spdmax_pos = [170 20 55 20];
  h_spdlbl_pos = [90 20 55 20];
  
  % add a button for play
  h_play = uicontrol(h_fig,...      
                     'callback','guislidr(gcf,''Play'');',...
                     'style','pushbutton',...
                     'string','Play',...
                     'position',h_play_pos);        

  % add a button for pause
  h_pause = uicontrol(h_fig,...     
                      'callback','guislidr(gcf,''Pause'');',...
                      'style','pushbutton',...
                      'string','Pause',...
                      'position',h_pause_pos);  
        
  % add a button for done
  h_done = uicontrol(h_fig,...      
                     'callback','guislidr(gcf,''Done'');',...
                     'style','pushbutton',...
                     'string','Done',...
                     'position',h_done_pos);        
        
  % add a button for stop (restart)
  h_stop = uicontrol(h_fig,...      
                     'callback','guislidr(gcf,''Restart'');',...
                     'style','pushbutton',...
                     'string','Restart',...
                     'position',h_stop_pos);        
        
  if version_number > 4
    h_sldr = uicontrol(h_fig,...
                     'callback','guislidr(gcf,''Slider Moved'');',...
                     'style','slider',...
                     'min',0,'max',2000,...
                     'value',300,...
                     'sliderstep',[0.005 .05],...        
                     'position',h_spdsld_pos);  
  else
    h_sldr = uicontrol(h_fig,...
                     'callback','guislidr(gcf,''Slider Moved'');',...
                     'style','slider',...
                     'min',0,'max',2000,...
                     'value',300,...
                     'position',h_spdsld_pos);  
  end % if version_number > 4                     

  speed_val = get(h_sldr,'value');

  h_min = uicontrol(h_fig,...
                    'style','text',...
                    'string','Slower',...
                    'BackgroundColor','w','ForegroundColor','k',...
                    'position',h_spdmin_pos);
  
  h_max = uicontrol(h_fig,...
                    'style','text',...
                    'string','Faster',... 
                    'BackgroundColor','w','ForegroundColor','k',...
                    'position',h_spdmax_pos);

  h_spd_lbl = uicontrol(h_fig,...
                    'style','edit',...
                    'callback','guislidr(gcf,''EditSpeed'');',...
                    'string',num2str(speed_val),... 
                    'BackgroundColor','w','ForegroundColor','k',...
                    'position',h_spdlbl_pos);

  play_val = get(h_play,'value');
  pause_val = get(h_pause,'value');
  done_val = get(h_done,'value');
  restart_val = get(h_stop,'value');
    
  handles = [h_sldr h_play h_pause speed_val play_val];
  handles = [handles pause_val done_val restart_val h_done h_stop h_spd_lbl];
  set(h_fig,'userdata',handles); 

% if the speed is being set externally
elseif strcmp(command_str,'SetSpeed')
  handles = get(gcf,'userdata');          % get the current userdata
  % set the speed value to the current value of the speed 
  % in the userdata (handles(2))
  set(handles(1),'value',handles(4));
  set(handles(11),'string',num2str(handles(4)));   

% if the speed is being set in the editable box in the GUI
elseif strcmp(command_str,'EditSpeed')
  % get the data from the user data property
  handles = get(gcf,'userdata');          % get the current userdata

  % get the current speed
  spd_string = get(handles(11),'string'); 

  % set the speed value to the current value of the speed 
  % in the userdata (handles(2))
  handles(4) = round(str2num(spd_string));
  set(handles(1),'value',handles(4));
  set(handles(11),'string',num2str(handles(4)));
  set(gcf,'userdata',handles);

% play button has been pushed
elseif strcmp(command_str,'Play')
  handles = get(gcf,'userdata');          % get the current userdata
  % set the play value to the current value of the play button (handles(2))
  handles(5) = get(handles(2),'value');   
  % set the pause value to 0
  handles(6) = 0;
  % reset all of the userdata with the new values
  set(gcf,'userdata',handles);

elseif strcmp(command_str,'Pause')
  handles = get(gcf,'userdata');     % get the current userdata
  % set the pause value to 1 (pause mode)
  handles(6) = 1;
  % set the play value to 0
  handles(5) = 0;
  % reset all of the userdata with the new values
  set(gcf,'userdata',handles);  

% the slider has been moved
elseif strcmp(command_str,'Slider Moved')
  handles = get(gcf,'userdata');     % get the current userdata
  % set the slider value to the current value of the slider bar (handles(1))
  handles(4) = round(get(handles(1),'value'));
  set(handles(1),'value',handles(4));
  set(handles(11),'string',num2str(handles(4)));   
  % reset all of the userdata with the new values
  set(gcf,'userdata',handles);

else

  handles = get(gcf,'userdata');
  h_sldr = handles(1);
  h_play = handles(2);
  h_pause = handles(3);
  speed_val = handles(4);
  done_val = handles(7); 
  restart_val = handles(8); 
  h_done = handles(9);
  h_stop = handles(10);
  h_spd_lbl = handles(11);
  
  speed_val = get(h_sldr,'value');
  play_val = get(h_play,'value');
  pause_val = get(h_pause,'value');
  done_val = get(h_done,'value');
  restart_val = get(h_stop,'value');
    
  handles = [h_sldr h_play h_pause speed_val play_val];
  handles = [handles pause_val done_val restart_val h_done h_stop h_spd_lbl];
  set(gcf,'userdata',handles);

end % if strcmp(command_str,'initialize')

%%%%% END ALGORITHM CODE %%%%%

% end of GUISLIDR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  MASK_STA function within ORB_ANIM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat, lon] = mask_sta(sta_loc,mask,alt,resolution);

% [lat, lon] = mask_sta(sta_loc,mask,alt,resolution);
%
% Computes azimuth/elevation mask projections to the surface of the Earth
% for ground stations. The resulting data is used in plotting ground station 
% masks on projections of the Earth.  
%   
% Input: 
%   sta_loc   - base station location in lat, lon, and hgt (1x3) [lat lon hgt]
%                lat and lon are in rad and height is in meters
%   mask - masking information (rad) (1x1), (nx3), or (nx4)
%        1x1 form is minimum elevation only [min_el]
%        nx3 form is min elevation and azimuth bounds [min_el start_az stop_az]
%        nx4 form is min and max elevation and azimuth bounds 
%          [min_el max_el start_az stop_az]
%        Azimuth start and stop are assumed to be clockwise.
%        Examples:
%         minimum elevation mask of 5 degrees (rad) (1x1) 
%              mask = .0873   
%         minimum elevation and azimuth bound triples (nx3)
%              mask = [.0873   pi/2   pi;       % 5 deg min el, 90->180 azimuth
%                      .1745   pi    pi/4]      % 10 deg min el, 180->90 azimuth
%                                               % wraps through 0
%   alt        - satellite altitude (m) (1x1)
%   resolution - azimuthal resolution of output points (rad) (1x1) (optional)
%                 default = .0175 (1 degree, n = 360)
% Output:
%   lat        - latitude for plotting ground station mask (nx1) (rad)
%   lon        - longitude for plotting ground station mask (nx1) (rad)
%
% See also VIS_DATA, ORB_ANIM 

% Written by: Jimmy LaMance 2/12/97
% Copyright (c) 1998 by Constell, Inc.

% Notes:
% this function will use the law of sines and cosines to compute ...
% 1) the slant range at the given elevation and satellite altitude
% 2) a set of n points (latitude and longitude) using azimuths from 
%    0:resolution:2*pi from the base station location

% WGS-84 constants
RADIUS_EARTH = 6378137.0;        % radius of the Earth WGS-84 value (m)

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
lon=[]; lat=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(3,4,nargin);
if ~isempty(msg)
  fprintf('%s  See help on MASK_STA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 4
  resolution = .0175;      % 1 deg
end

estruct.func_name = 'MASK_STA';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'sta_loc';
estruct.variable(1).req_dim = [1 3];
estruct.variable(1).var = sta_loc;
  
estruct.variable(2).name = 'mask';
estruct.variable(2).req_dim = [1 1; 901 3; 901 4];
estruct.variable(2).var = mask;
  
estruct.variable(3).name = 'alt';
estruct.variable(3).req_dim = [1 1];
estruct.variable(3).var = alt;
  
estruct.variable(4).name = 'resolution';
estruct.variable(4).req_dim = [1 1];
estruct.variable(4).var = resolution;

% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% Set up value for wraps in longitude.  The default will be wrapping at
% 0 and 360 degrees.  At these points, a new line of output will be generated.
long_wrap = [0 2*pi];

% Fill in the mask variable to be nx4
if size(mask,2) == 1
  mask = [mask(:,1) pi/2*ones(size(mask,1),1) ...
          zeros(size(mask,1),1) 2*pi*ones(size(mask,1),1)];
elseif size(mask,2) == 3
  mask = [mask(:,1) pi/2*ones(size(mask,1),1) mask(:,2:3)];
end

% compute base station coordinates in ECEF
base_ecef = lla2ecef(sta_loc);           

% compute base station height from CM of the Earth
[base_norm, base_h] = normvect(base_ecef);

% initialize lat and lon output matrices
all_data.lon = [];
all_data.lat = [];

% 1 degree tolerance in azimuth for drawing masks, anything less than 1 degree
% will be not be drawn
az_tol = 1 * pi/180;       

% Define 2 type of masks, complete circles and partial circles
% each type will be draw with the same routine with the excpetion of
% drawing the edge lines of constant azimuth
I_full_circ = find(rem(abs(mask(:,4) - mask(:,3)),2*pi) < az_tol);
I_part_circ = find(rem(abs(mask(:,4) - mask(:,3)),2*pi) > az_tol);

% save all of the masking information
mask_all = mask;

mask_counter = 1;

% start with the partial circles
if any(I_part_circ)
  mask = mask_all(I_part_circ,:);

  % sort the coverage masks such that the starting azimuth is always increasing
  if size(mask,1) > 1
    [sorted_az, I_sort] = sort(mask(:,3));
    mask = mask(I_sort,:);
  end % if size(mask,1) > 1

  % Initialize the variables for temporary storage of the lat/lon data
  % until a full set to generate a patch is filled.
  lon_temp = [];
  lat_temp = [];                
  first_flag = 1;   % used to draw the first side of the mask
  
  % loop over the sets of azimuth/elevation quadrulpes
  for i = 1:size(mask,1)

    % Find out which mask is the next
    if i == size(mask,1)
      next_mask = 1;
    else
      next_mask = i + 1;
    end % if i == size(mask,1)
    
    % Check to see if this mask (mask(i,:)) is continuous with the i+1 mask.
    % If the masks are not next to each other, a dummy mask will be filled
    % in with a minimum elevation of pi/2.
    if abs(mask(next_mask,3) - mask(i,4)) > az_tol
      continuous_flag = 0;
    else
      continuous_flag = 1;
    end % if abs(mask(next_mask,3) - mask(i,4)) > az_tol
    
    % Start drawing lines....
    % If the first flag is set, then we have to draw a constant azimuth
    % line on the min azimuth side of the mask from the max to min elevations
    if first_flag == 1;                                                      
      % Check to see if the last mask is continutous with this one
      last_mask = size(mask,1);
      if abs(mask(last_mask,4) - mask(i,3)) > az_tol   
        elevation = [mask(i,2):-resolution:mask(i,1)]';
      else             
        start_mask = min([mask(i,1) mask(last_mask,1)]);
        stop_mask = max([mask(i,1) mask(last_mask,1)]);
        elevation = [start_mask:resolution:stop_mask]';
      end
      azimuth = mask(i,3);

      [lat_this_line, lon_this_line] = ...
              az_line(elevation,azimuth,sta_loc,base_h,alt,resolution);
        
      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];      
      first_flag = 0;
    end

    % The next line is the min elevation line across the azimuth range
    elevation = mask(i,1);
    if mask(i,3) < mask(i,4)                     % no azimuth wrapping
      azimuth = [mask(i,3):resolution:mask(i,4)]';
    else                                         % azimuth wrapping
      azimuth = [mask(i,3):resolution:mask(i,4)+2*pi]';
      azimuth = rem(azimuth,2*pi);                       
    end % if mask(i,3) < mask(i,4)               
 
    [lat_this_line, lon_this_line] = ...
            el_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

    lon_temp = [lon_temp; lon_this_line];
    lat_temp = [lat_temp; lat_this_line];

    % If the continuous flag is set, then we have to draw a constant azimuth
    % lines on the max azimuth side of the mask to connect with the mext 
    % set of masking information
    if continuous_flag == 1;                          
      if mask(i,1) < mask(next_mask,1)                            
        elevation = [mask(i,1):resolution:mask(next_mask,1)]';
      else
        elevation = [mask(next_mask,1):resolution:mask(i,1)]';
      end % if mask(i,1) < mask(next_mask,1)            

      azimuth = mask(i,4);
      
      [lat_this_line, lon_this_line] = ...
              az_line(elevation,azimuth,sta_loc,base_h,alt,resolution);
        
      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];
      
      % Now for the second line connecting the max elevations
      if mask(i,2) < mask(next_mask,2)                            
        elevation = [mask(i,2):resolution:mask(next_mask,2)]';
      else
        elevation = [mask(next_mask,2):resolution:mask(i,2)]';
      end % if mask(i,2) < mask(next_mask,2)            
        
      azimuth = mask(i,4);

      [lat_this_line, lon_this_line] = ...
              az_line(elevation,azimuth,sta_loc,base_h,alt,resolution);
        
      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];

    % If the continous flag is not set, we draw a single line
    else   
      elevation = [mask(i,1):resolution:mask(i,2)]';
      azimuth = mask(i,4);
      
      [lat_this_line, lon_this_line] = ...
              az_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];
    end

    % The next line is the max elevation line across the azimuth range
    elevation = mask(i,2);
    if abs(elevation - pi/2) > az_tol
      if mask(i,3) < mask(i,4)                     % no azimuth wrapping
        azimuth = [mask(i,4):-resolution:mask(i,3)]';
      else                                         % azimuth wrapping
        azimuth = [mask(i,4):-resolution:mask(i,3)+2*pi]';
        azimuth = rem(azimuth,2*pi);                       
      end % if mask(i,3) < mask(i,4)               
 
      [lat_this_line, lon_this_line] = ...
              el_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];
    end % if abs(elevation - pi/2) > az_tol
    
    % If this mask is continuous with the next one, don't increment
    % the mask_counter to create a new structure element of lat/lon data.  If it
    % is not continuous, increment the counter, store the lat and lon
    % temp values into the current counter, and reset the
    % temporary values to empty matrices.

    if continuous_flag == 0 |  i == size(mask,1)
      all_data(mask_counter).lon = lon_temp;  
      all_data(mask_counter).lat = lat_temp;  
      mask_counter = mask_counter + 1; 
      first_flag = 1;
      lon_temp = [];
      lat_temp = [];
    end
      
  end % for i = 1:size(mask,1)
end % if any(I_part_circ)

% now repeat the process for circular masks
if any(I_full_circ)
  mask = mask_all(I_full_circ,:);

  % Initialize the variables for temporary storage of the lat/lon data
  % until a full set to generate a patch is filled.
  lon_temp = [];
  lat_temp = [];                

  % loop over the sets of azimuth/elevation triples
  for i = 1:size(mask,1)        

    if mask(i,3) < mask(i,4)                     % no azimuth wrapping
      azimuth = [mask(i,4):-resolution:mask(i,3)]';
    else                                         % azimuth wrapping
      azimuth = [mask(i,4):-resolution:mask(i,3)+2*pi]';
      azimuth = rem(azimuth,2*pi);                       
    end % if mask(i,3) < mask(i,4)               
 
    % The first line is the min elevation line across the azimuth range
    elevation = mask(i,1);
    if abs(elevation - pi/2) > az_tol
      [lat_this_line, lon_this_line] = ...
              el_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];
    end % if abs(elevation - pi/2) > az_tol

    % The next line is the max elevation line across the azimuth range
    elevation = mask(i,2);
    if abs(elevation - pi/2) > az_tol
      [lat_this_line, lon_this_line] = ...
              el_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

      lon_temp = [lon_temp; lon_this_line];
      lat_temp = [lat_temp; lat_this_line];
    end % if abs(elevation - pi/2) > az_tol

    all_data(mask_counter).lon = lon_temp;  
    all_data(mask_counter).lat = lat_temp;  
    mask_counter = mask_counter + 1; 
    first_flag = 1;
    lon_temp = [];
    lat_temp = [];

  end % for i = 1:size(mask,1)
end % if any(I_full_circ)

% Get the structure data into an array format.  Start by finding the
% maximum size of the entire structure array.  At this time the only
% way to do this is to check each structure element.  This is not to bad
% as far as vectorization because there should only be a few
% structure elements to loop over.
max_size = 0;
for i = 1:size(all_data,2)
  max_size = max([max_size length(all_data(i).lon)]);
end

% Set up an empty matrix to the max size
lon_new = ones(max_size,size(all_data,2)) * NaN;
lat_new = ones(max_size,size(all_data,2)) * NaN;

% Now fill in the new lat and lon matrices
for i = 1:size(all_data,2)
  lon_new(1:length(all_data(i).lon),i) = all_data(i).lon;
  lat_new(1:length(all_data(i).lat),i) = all_data(i).lat;
end
                      
% Clear out the old lat and lon
lon = lon_new;
lat = lat_new;
                      
I = find(lon < 0);
if any(I)
  lon(I) = lon(I) + 2 * pi;
end % if any(I)

return
keyboard
% Now that we have all of the lat and lon points, break them up into rows
% where there is no wrap around the longitude wrapping points.  
I_wrap = find(abs(diff(lon)) > 3*pi/2);

if ~isempty(I_wrap)                        
  I_wrap = [I_wrap; size(lon,1)];
  new_lon = ones(size(lon,1),length(I_wrap)) * NaN;
  new_lat = ones(size(lon,1),length(I_wrap)) * NaN;
  start_index = 1;
  
  for i = 1:length(I_wrap)
    new_lon(start_index:I_wrap(i),i) = lon(start_index:I_wrap(i));
    new_lat(start_index:I_wrap(i),i) = lat(start_index:I_wrap(i));
    start_index = start_index + I_wrap(i);
  end

  lat = new_lat;
  lon = new_lon;
end

%%%%% END ALGORITHM CODE %%%%%

% end of MASK_STA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine within MASK_STA
function [lat_temp, lon_temp] = ...
       az_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

% WGS-84 constants
RADIUS_EARTH = 6378137.0;        % radius of the Earth WGS-84 value (m)

% compute the angle (B) between the slant range vector and the 
% satellite range vector using the law of sines                                              
B = asin(base_h * sin(pi / 2 + elevation) / (RADIUS_EARTH + alt));

% compute the interior angle (A) at CM of the Earth
A = pi - B - (pi/2 + elevation);

% compute the slant range using the law of cosines
slant = sqrt((RADIUS_EARTH + alt)^2 + base_h^2 - ...
              2 * (RADIUS_EARTH + alt) * (base_h) * cos(A));

range(:,1) = slant .* cos(elevation) .* cos(azimuth);
range(:,2) = slant .* cos(elevation) .* sin(azimuth); 
range(:,3) = -slant .* sin(elevation); 

% Convert the range vector to ECEF
range_ecef = ned2ecef(range,sta_loc);      

base_ecef = lla2ecef(sta_loc);

% Add the range to the base station to get ECEF positions to be projected
% on the surface of the Earth.
sat_pos_ecef = ones(size(range_ecef,1),1) * base_ecef + range_ecef;

% Convert to lla coordinates
sat_pos_lla = ecef2lla(sat_pos_ecef); 

lat_temp = sat_pos_lla(:,1);
lon_temp = sat_pos_lla(:,2);

% end of function az_line       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine within MASK_STA
function [lat_temp, lon_temp] = ...
       el_line(elevation,azimuth,sta_loc,base_h,alt,resolution);

% WGS-84 constants
RADIUS_EARTH = 6378137.0;        % radius of the Earth WGS-84 value (m)

% compute the angle (B) between the slant range vector and the 
% satellite range vector using the law of sines                                              
B = asin(base_h * sin(pi / 2 + elevation) / (RADIUS_EARTH + alt));

% compute the interior angle (A) at CM of the Earth
A = pi - B - (pi/2 + elevation);

% compute the slant range using the law of cosines
slant = sqrt((RADIUS_EARTH + alt)^2 + base_h^2 - ...
              2 * (RADIUS_EARTH + alt) * (base_h) * cos(A));

% compute the range vectors for these elevation changes    
range(:,1) = slant .* cos(elevation) .* cos(azimuth);
range(:,2) = slant .* cos(elevation) .* sin(azimuth); 
range(:,3) = -slant .* sin(elevation); 

% Convert the range vector to ECEF
range_ecef = ned2ecef(range,sta_loc);      

base_ecef = lla2ecef(sta_loc);

% Add the range to the base station to get ECEF positions to be projected
% on the surface of the Earth.
sat_pos_ecef = ones(size(range_ecef,1),1) * base_ecef + range_ecef;

% Convert to lla coordinates
sat_pos_lla = ecef2lla(sat_pos_ecef); 

lat_temp = sat_pos_lla(:,1);
lon_temp = sat_pos_lla(:,2);

% end of function az_line