function [M, color_map] = ...
          ex_dop_m(start_utc,duration,time_step,spatial_res,dop_type,mask,alm_file);

% [M, color_map] = ...
%         ex_dop_m(start_utc,duration,time_step,spatial_res,dop_type,mask,alm_file);
%                          
% Example function to generate a movie of GPS DOP values.  This movie shows
% the user chosen DOP as a function of time over the surface of the Earth. 
% The output is a movie matrix (M) and the associated color map such that 
% the movie can be saved using the WRITEMOV or SAVE functions or it can be
% replayed using the PLAYMOV function.
%
% Input:  
%   start_utc   - UTC start time [year month day hour minute second] (1x6)
%                  (optional), if not provided, the computer time is used
%   duration    - total time duration (1x1) (sec) (optional) default = 3600
%   time_step   - step size for computation of visibility (1x1) (sec) 
%                  (optional), default = 600. Smaller time steps take more 
%                  memory, but do not significantly increase the run time.
%   spatial_res - spatial resolution of the DOP computations.  (1x1) (deg)
%                  (optional), default = 5.  Smaller values of the spatial
%                  resolution provide smoother movies, but take longer to 
%                  generate.
%   dop_type    - type of DOP to show in the movie (optional) (1x1)
%                  1 - GDOP, 2 - PDOP, 3 - HDOP, 4 - VDOP, 5 - TDOP
%                  default = 1, GDOP
%   mask        - elevation mask_vis (rad) (1x1 or nx3) (optional) default = 0
%                  see help on VIS_DATA for more detail about mask parameters
%   alm_file    - almanac file name to read (nx1) (string) (optional)
%                  default will find the GPS almanac most recent to the start time 
% Output: 
%   M           - Matlab movie matrix, size is dependent on the pixel resolution
%                  of the figure window and the number of time steps 
%   color_map   - color map corresponding to the movie matrix
%
% Note: For small spatial resolutions, it may take a significant period of time
%       to generate the output.  In these cases, the screen saver should be 
%       disabled or the resulting movie matrix will only contain screen saver
%       images.  The movie that is generated is a screen capture of the
%       active movie window space.  See help on MOVIE for more details.
%
% See also PLAYMOV, WRITEMOV, MOVIE

% Written by: Jimmy LaMance 2/19/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: UTC2GPS, LLA2ECEF, READYUMA, ALM2GEPH, PROPGEPH,
%                   LOS, ECEF2NED, NED2AZEL, VIS_DATA, NED2DOPS,
%                   GPS2UTC, PLAYMOV

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% The level of error checking in the example functions is less stringent than
% the standard GNSS Toolbox functions. 

% check inputs
if nargin < 1
  % get the data from the computer clock
  start_utc = clock;  
end % if nargin < 1

% if the duration is not given, set to the default value
if nargin < 2
  duration = 3600;   % set the duration to the default
end % if nargin < 2                                    

% if the time step is not given, set to the default value
if nargin < 3
  time_step = 600;   % set the time step to the default
end % if nargin < 3  

% if the spatial resolution is not given, set to the default value
if nargin < 4 
  spatial_res = 15;   % set the default spatial resolution
end % if nargin < 4

% if the DOP type is not given, set to the default value
if nargin < 5
  dop_type = 1;   % set the default DOP type, GDOP
end % if nargin < 4

% if there is no mask input, set to the default value
if nargin < 6
  mask = 0;
end % if nargin < 6                          

% Set the almanac file to read to be the most recent GPS almanac
if nargin < 7
  alm_file = find_alm;
end % if nargin < 7  

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

colordef white
% set a maximum DOP value to show on the plot
maximum_dop_to_show = 20;

% compute the stop time based on duration
% first convert the start time to GPS time
start_gps = utc2gps(start_utc);

% add on the duration
stop_gps(1) = start_gps(1);
stop_gps(2) = start_gps(2) + duration;

% check for week roll-overs
if stop_gps(2) > 604800
  stop_gps(1) = stop_gps(1) + fix(stop_gps(2) / 604800);
  stop_gps(2) = start_gps(2) + rem(stop_gps(2), 604800);
end % if stop_gps(2) > 604800

% conversion from degrees to radians
d2r = pi / 180;

% generate the lat/lon/time points for the positions on the Earth
lats = [-90:spatial_res/2:90]';
lons = [0:spatial_res:360]';
num_lats = size(lats,1);
num_lons = size(lons,1);
num_points = num_lats * num_lons;

% set the height field for all locations to zero
hgt = zeros(num_points,1);

% generate indices into the lat/lons for setting up a single computation of 
% fixed positions on the Earth.  The way this will work is as follows...
% 1) the index to the resulting points will be 
%       1:num_points = [1 2 3 ... num_points]
% 2) generate an index into the latitudes for each of the num_points as
%       [1 2 3 ...num_lons 1 2 3 ... num_lons ... num_lons]
%       where the sequence repeats for the number of lons
% 3) generate an index into the longitudes for the of the num_points as
%       [1 1 1 1 1... 1 2 2 2 2 2... 2 3 3 3 3...num_lats]
%       where there are as many 1s, 2s, 3s, etc as there are lats

% step #1, index to the output points
point_index = [1:num_points]';          

% step #2, index to the latitudes
lat_index = zeros(size(point_index,1),1);

for i = 1:num_lons
  lat_index((i-1)*num_lats+1:i*num_lats) = 1:num_lats;
end

% step #3, index to the longitudes
% start by putting 1s at each new latitude point (ie where the change from
% 1 to 2, 2 to 3, etc occurs and filling the remainder with zeros
lon_index = zeros(size(point_index,1),1);                        
lon_index(1:num_lats:num_points) = ones(size([1:num_lats:num_points]',1),1);

% fill in the rest of the matrix with the appropriate values using cumsum
lon_index = cumsum(lon_index);

% build the ECEF position vector for the lat/lon values using the
% indices computed above (this is a vectorized method)
x_pos_all = lla2ecef([lats(lat_index) * d2r, lons(lon_index) * d2r, hgt]);

% Read in the almanac file
[gps_alm] = readyuma(alm_file);

% sort out the unhealthy satellites
I_gps_good = find(gps_alm(:,2) == 0);
gps_alm = gps_alm(I_gps_good,:);

% convert the almanacs to ephemeris format
[gps_ephem] = alm2geph(gps_alm);

% how many satellites are in the alamanc/ephemeris
total_gps_satellites = size(gps_ephem,1);

% build up a matrix of times to line up with the Earth based positions
% computed above
t_user(:,2) = [start_gps(2):time_step:stop_gps(2)]';
t_user(:,1) = ones(size(t_user,1),1) * start_gps(1);

% compute the number of times
num_time_intervals = size(t_user,1);

% adjust the t_user times for week rollovers
I = find(t_user(:,2) > 604800);
if any(I)
  t_user(I,1) = rem(t_user(I,2), 604800) + start_gps(1);
  t_user(I,2) = t_user(I,2) - (t_user(I,1) - start_gps(1)) * 604800;
end % if

% compute satellite positions in ECEF at the times specified by the inputs
[t_gps, prn_gps, x_gps, v_gps] = propgeph(gps_ephem, t_user);

% set a timer, a later call to the clock function will time how long it 
% takes for the DOPs computation
t0 = clock;                                            

% allocate the DOP matrix
dop = zeros(num_points,num_time_intervals);

% loop over all the points in the Earth to be evaluated
% the time sequence at each location will be computed in the loop
% so there is no outer time loop
for i = 1:num_points

  [t_los_gps, gps_los] = los(t_gps(1,:), x_pos_all(i,:), t_gps, [prn_gps x_gps]);

  % convert LOS in ECEF to North-East-Down (NED)
  [gps_los_ned] = ecef2ned(gps_los, ecef2lla(x_pos_all(i,:)));

  % Compute az and el
  [az el] = ned2azel(gps_los_ned);
  
  % Compute masking
  [az_el, I_pass] = vis_data(mask,[az el]);
  
  % compute DOPS
  [gps_dops, t_gps_dops, num_gps_sats] = ...
         ned2dops(gps_los_ned(I_pass,:),t_los_gps(I_pass,:));

  % verify that the dops are filled correctly. NED2DOPS does not return
  % any dops values when there are no satellites visibles.  It's possible
  % that the gps_dop matrix may be empty or partially filled.
  if size(gps_dops,1) ~= num_time_intervals

    % create a temporary variable to store data
    temp_dops = ones(num_time_intervals,5) * inf;
    
    % check to see if there were any valid dops returned from NED2DOPS
    if ~isempty(gps_dops)
     
      % loop over the number of times and fill in the DOPs if required
      for j = 1:num_time_intervals        
        I_this_time = find(t_gps_dops(:,2) == t_user(j,2));
       
        % if there is output data for this time, fill in the temp matrix with it
        % otherwise, the inf value that temp_dops was initialized to will 
        % be maintained
        if ~isempty(I_this_time)
          temp_dops(j,:) = gps_dops(I_this_time,:);
        end % if ~isempty(I_this_time)
      
      end % for j = 1:num_time_intervals           
    
      % now that the temporary variable is filled correctly, rename it
      gps_dops = temp_dops;
    
      clear temp_dops I_this_time  
    
    else   % NED2DOPS returned all blanks
      
      % set the output to the temp values of inf
      gps_dops = temp_dops;
                                                
    end % if ~isempty(t_gps_dops)
                                                  
  end % if size(gps_dops,1) ~= num_time_intervals

  % fill in the dop matrix that will be used to generate the output plots  
  dop(i,:) = gps_dops(:,dop_type)';                                      
  
end % for  

% clear out variables that are no longer needed
clear t_los_gps prn_gps_los gps_los gps_dops t_gps_dops num_gps_sats  
clear x_user_all_gps
clear t_los_gps gps_los x_user_all_gps prn_gps_los x_pos t_gps

% compute the elapsed time in the DOPs routine and print that 
% information to the screen
t1 = etime(clock,t0);
fprintf('Elapsed time was %6.2f minutes computing DOPs\n',t1 / 60);
fprintf('for %d points on the Earth and %d time steps for a\n',...
         num_points,num_time_intervals);                   
fprintf('total of %d DOP computations.\n',num_points * num_time_intervals);

% begin the work for generating the figures that will comprise the movie

% load the topo.mat file (comes with Matlab) that will be used to generate
% the land outlines, the 0 contour line of the topo is the boundary between
% the ocean and the land 
load topo

% make a movie of the DOPS   
% first generate the axis data and
% scale the DOP to make it fit a 64 color map
min_dop = min(min(dop));
max_dop = max(max(dop));

min_value = floor(min_dop);
max_value = ceil(min([max_dop,maximum_dop_to_show]));

% put ticks on the color bar at interger values
num_ticks = max_value - min_value + 1;

% scale all the values between 0 and 64
to = fix(64*(dop-min_value)/(max_value-min_value)+1);

% reshape the scaled dops matrix for input to the image function
dop_now = reshape(to(:,1),num_lats,num_lons);

% check the current version of Matlab in use
matlab_version = version;      % full version string

% get the first digit (either 4 or 5), or report an error
version_number = str2num(matlab_version(1));

% check that this is version 4 or 5, if not report an unsupported
% version and return
if version_number ~= 4 & version_number ~= 5
  fprintf('This function is only supported under Matlab versions 4 and 5.\n') ;
  return
end % if version_number ~= 4 & version_number ~= 5
                                 
% generate a new figure
figure

% set the figure background color to white
set(gcf,'color','w');

% put the DOP values on the figure using the image function
img_handle = image([1 360],[-90 90],dop_now);, axis xy, axis image

% turn the hold on so the contour lines will be overlaid on the image
hold on

% generate the contour plot
% check the Matlab version because the contour command changed between
% versions 4 and 5
if version_number < 5
  % version 4 command
  [c, topo_handles] = contour(1:360,-89:90,topo,[0 0]); 
else
  % version 5 command
  [c, topo_handles] = contour(1:360,-89:90,topo,[0 0],'k-'); 
end % if version_number < 5

% set the property for the contour lines so that they will not be earsed
set(topo_handles,'EraseMode','none'); 

% set the axis for display of the world
set(gca,'XLim',[0 360],'YLim',[-90 90]);

% make the lines black and specify a width
set(topo_handles,'Color','k','LineWidth',1.25);

% add the color bar and set some of it's attributes
cb_handle = colorbar('horiz');       % make the bar horizontal
set(cb_handle,'XTickMode','manual','XTickLabelMode','manual');
set(cb_handle,'YTickMode','manual','YTickLabelMode','manual');

% compute the color resolution between integer DOP values
xt_resolution = 65 / (num_ticks-1);

% add 1 to the integer count to put the end tick mark on the colorbar
xt = [1:xt_resolution:64];

% make sure that the last value in xt is 64
if xt(length(xt)) ~= 64
  xt = [xt 64];
end % if xt(length(xt)) ~= 64

% set the limits on the color bar to be the same as the steps
% this is done because the last color may be slightly less than 64
set(cb_handle,'XLim',[xt(1) xt(length(xt))]);

% set the axis tick labeling mode to manual
set(gca,'XTickMode','manual','XTickLabelMode','manual');
set(gca,'YTickMode','manual','YTickLabelMode','manual');

% define the X and Y tick locations and labels for the world plot
str1 = str2mat('0','30','60','90','120','150','180','210');
str2 = str2mat('240','270','300','330','360');
x_tick_marks = str2mat(str1,str2);
x_tick_vals = [0:30:360]';                        

y_tick_marks = str2mat('-90','-60','-30','0','30','60','90');
y_tick_vals = [-90:30:90]';                        

% define the tick labels for the color bar
% make the last color ba tick label be greater than (>) the last
% value.  this is done because there is a hard limit of a DOP of 20
% in the plotting.  To change the hard stop at 20, modify the variable
% maximum_dop_to_show at the beginning of the function. 
if version_number < 5
  % version 4 command
  x_tick_labels(1:num_ticks-1,:) = str2mat(num2str([min_value:max_value-1]))'; 
else
  % version 5 command
  x_tick_labels(1:num_ticks-1,:) = str2mat(num2str([min_value:max_value-1]')); 
end % if version_number < 5

last_tick = sprintf('>%d',max_value);
x_tick_labels = str2mat(x_tick_labels, last_tick);

% set the tick labels based on which version of Matlab is in use
if version_number == 4   
  set(cb_handle,'XTick',xt,'XTickLabels',x_tick_labels);
  set(gca,'XTick',x_tick_vals,'XTickLabels',x_tick_marks);
  set(gca,'YTick',y_tick_vals,'YTickLabels',y_tick_marks);
else   % version 5
  set(cb_handle,'XTick',xt,'XTickLabel',x_tick_labels);
  set(gca,'XTick',x_tick_vals,'XTickLabel',x_tick_marks);
  set(gca,'YTick',y_tick_vals,'YTickLabel',y_tick_marks);
end % if version_number == 4 

% add x and y labels to the world plot
xl = xlabel('Longitude','erasemode','none');
yl = ylabel('Latitude','erasemode','none');

% add title with the date/time stamp
% start with a string to identify which DOP is being plotted,
% the order corresponds to the dop_type variable (dop_string(1,:) = 'GDOP')
dop_string = str2mat('GDOP','PDOP','HDOP','VDOP','TDOP'); 

% label the color bar with the appropriate DOP value
cb_xstring = sprintf('%s Value',dop_string(dop_type,:));
cb_x_label_handle = get(cb_handle,'XLabel');
set(cb_x_label_handle,'String',cb_xstring);

% compute UTC from GPS time
utc_time = gps2utc(t_user(1,:));    

% set the date string
date_text = sprintf('%d/%d/%d  %d:%02d:%02.2f',...   
            utc_time(2),utc_time(3),utc_time(1)-1900,utc_time(4:6));
            
% incorporate the date text into the title string            
title_text = sprintf('%s for %s',dop_string(dop_type,:),date_text);

% add the title to the figure
title(title_text)

% make this figure fill the entire screen
% set the base units to pixels
set(0,'Units','pixels');    

% get the screen size from the base properties
screen_size = get(0,'ScreenSize');   

% amount to adjust the screen size by to get the window to fit (pixels)
adjust = [60 60 -120 -120];         

% set the figure position to fill the screen
set(gcf,'Position',screen_size + adjust);    

% allocate the movie matrix, this greatly speeds up the movie generation
% because the M matrix can become very large
M = moviein(num_time_intervals,gcf);

%keyboard
% get the current figure as the first frame of the movie
M(:,1) = getframe(gcf);

% turn the hold off so that the figure will be replaced by the next frame
hold off

% now generate the remainder of the frames
for i = 1:num_time_intervals

  % get the DOPs into the right format and generate the image
  dop_now = reshape(to(:,i),num_lats,num_lons);
  img_handle = image([1 360],[-90 90],dop_now);, axis xy, axis image

  % turn to hold on so we can overlay the topography contours
  hold on
  
  % generate the contour plot
  if version_number < 5
    [c, topo_handles] = contour(-179:180,-89:90,topo,[0 0]); 
  else
     [c, topo_handles] = contour(1:360,-89:90,topo,[0 0],'k-'); 
  end % if version_number < 5

  set(topo_handles,'EraseMode','none'); 

  % set the axis for display of the world
  set(gca,'XLim',[1 360],'YLim',[-90 90]);

  % make the topo lines black and specify a width
  set(topo_handles,'Color','k','LineWidth',1.25);

  % add the color bar
  cb_handle = colorbar('horiz');       % make the bar horizontal
  set(cb_handle,'XLim',[xt(1) xt(length(xt))]);
  set(cb_handle,'XTickMode','manual','XTickLabelMode','manual');

  if version_number == 4   
    set(cb_handle,'XTick',xt,'XTickLabels',x_tick_labels);
    set(gca,'XTick',x_tick_vals,'XTickLabels',x_tick_marks);
    set(gca,'YTick',y_tick_vals,'YTickLabels',y_tick_marks);
  else   % version 5
    set(cb_handle,'XTick',xt,'XTickLabel',x_tick_labels);
    set(gca,'XTick',x_tick_vals,'XTickLabel',x_tick_marks);
    set(gca,'YTick',y_tick_vals,'YTickLabel',y_tick_marks);
  end % if version_number == 4 
  
  % label the x-axis of the color bar
  cb_x_label_handle = get(cb_handle,'XLabel');
  set(cb_x_label_handle,'String',cb_xstring);
  
  % add x and y labels
  xl = xlabel('Longitude','erasemode','none');
  yl = ylabel('Latitude','erasemode','none');

  % add title with the date/time stamp
  % compute UTC from GPS time
  utc_time = gps2utc(t_user(i,:));     
  
  % set the date string
  date_text = sprintf('%d/%d/%d  %d:%02d:%02.2f',...  
            utc_time(2),utc_time(3),utc_time(1)-1900,utc_time(4:6));
  
  title_text = sprintf('%s for %s',dop_string(dop_type,:),date_text);
  title(title_text)
  
  % get this figure as a frame in the movie
  M(:,i) = getframe(gcf);

  % turn hold off, so that the next frame replace this one  
  hold off
end

% get the color map used for the movie so the replay will use the same colors
color_map = colormap;

% close the figure used to generate all of the plots
close(gcf)

% play the movie by opening a figure with the correct dimensions using the
% PLAYMOV command.  This command generates a figure with the correct dimensions,
% play the movie 4 times at 2 Hz, with the correct color map
playmov(M,4,2,color_map);
 
% end of EX_DOP_M
