function [fig_handles] = makeplot(visible_data, pass_numbers, num_sats_vis, ...
                                  gps_dops, label_string, ...
                                  plot_selection, min_elev, vertical_offset)

% fig_handles = makeplot(visible_data, pass_numbers, num_sats_vis, ...
%                                  gps_dops, label_string, ...
%                                  plot_selection, min_elev, vertical_offset);
%
% Function to create the five most used plots, azimuth, elevation, sky plot,
% number of visible satellites, and DOPS for a single object.
%
% Inputs:
%   visible_data    - [azimuth, elevation, GPS week, GPS seconds, PRN] that are
%                     visible (nx5) (radians)
%   pass_numbers    - pass number of each observation (nx1)
%   num_sats_vis    - [GPS week, GPS Sec, number of visible satellites] (mx3)
%   gps_dops        - [GPS week, GPS Sec, GDOP, PDOP, HDOP, VDOP, TDOP] (jx7)
%   label_string    - string identifying object to be plotted (1x??)
%                     i.e. 'Boulder' for ground site at Boulder
%                       'User Satellite 1', etc.
%   plot_selection - (1x9) optional array to limit the plots that are created
%                    each row [1=plot sky plot, 1=plot elevation, 1=plot azimuth, 
%                              1=plot # of visible satellites, 1=plot PDOP,
%                              1=plot GDOP, 1=plot HDOP, 1=plot VDOP, 1=plot TDOP]
%                     (optional), default = all plots
%   min_elev        - Minimum elevation to be plotted (radians) (optional)
%                      default = -pi/2 (all data plotted)
%   vertical_offset - objects are plotted above the previously plotted object. (1x1)
%                     this integer identifies the number of plots to shift upward.
%                     (optional), default = 0.
%
% See also PLOTPASS, PLOTSKY

% Written by: Maria J. Evans May 98
% Copyright (c) 1998 by Constell, Inc.

% functions called: PLOTPASS, GPST2SEC, GPS2UTC, PLOTSKY

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize output variables
fig_handles=[];

if nargin < 6
  plot_selection = ones(1,9);
end

if nargin < 7
  min_elev = min(visible_data(:,2));
end

if nargin < 8 
  vertical_offset = 0;
end

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% Plotting from here on.
    
% colordef none

% Get the screen size in pixels to use for location of plots
set(0,'units','pixels');
screen_size = get(0,'screensize');
y_max = screen_size(2) + screen_size(4) - 60;
x_max = screen_size(1) + screen_size(3) - 70;
x_step = 80;
y_step = 60;

  % Determine the location for the first azimuth plot in upper right corner
  x_min = x_max/2;
  y_min = min(30+y_step*(vertical_offset), 30+50*(vertical_offset));

  % For each location, compute the plots
  n_plots = 0;
  DOP_label = ['GDOP';'PDOP';'HDOP';'VDOP';'TDOP';];

  if plot_selection(1,3)==1,
    % create a figure for the azimuth plot
    n_plots = n_plots + 1;

    % create a title for the azimuth plot
    title_str_az = sprintf('GPS Azimuths for %s', label_string);

    plot_data(:,n_plots) = visible_data(:,1)*180/pi;
    plot_title(n_plots,1:length(title_str_az)) = title_str_az;
    y_label(n_plots,1:13) = 'Azimuth (deg)';
    y_scale(n_plots,:) = [0 360];

  end;

  if plot_selection(1,2)==1,

    % create a figure for the elevation plot and shift it left of the azimuth plot
    n_plots = n_plots + 1;

    % create a title for the elevation plot
    title_str_el = sprintf('Elevations for %s', label_string);

    plot_data(:,n_plots) = visible_data(:,2)*180/pi;
    plot_title(n_plots,1:length(title_str_el)) = title_str_el;
    y_label(n_plots,1:15) = 'Elevation (deg)';
    y_scale(n_plots,:) = [min_elev*180/pi 90];
  end;

  if n_plots > 0,  % Call plotpass to plot az and el
    % generate the plot of elevation versus time
    fig_handles(1:n_plots) = plotpass(visible_data(:,3:4),plot_data, ...
                  [visible_data(:,5), pass_numbers], ...
                  plot_title,y_label,y_scale);
    set(fig_handles(1),'position',[x_min y_min x_max/2 y_max/2]);
    if n_plots==2,
      set(fig_handles(2),'position',[x_min-x_step y_min x_max/2 y_max/2]);
    end;
    clear plot_data
  end;

  if plot_selection(1,4)==1,

    % create a figure for the number of visible satellites plot and shift it left of the elevation plot
    n_plots = n_plots + 1;

    % create a title for the visible satellites plot
    title_str_vis = sprintf('Number of Visible Satellites for %s', label_string);

    % generate the data for plot of visible satellites versus time
    prn_vis = zeros(size(num_sats_vis,1),1);
    y_label=['# of Visible Satellites'];
    fig_handles(n_plots) = plotpass(num_sats_vis(:,1:2), num_sats_vis(:,3), prn_vis, ...
                       title_str_vis, y_label, ...
                       [min(num_sats_vis(:,3))-1, max(num_sats_vis(:,3))+1]);
    set(fig_handles(n_plots),'position',[x_min-x_step*2 y_min x_max/2 y_max/2]);
  end;

  % Check how many DOPs are plotted for this observer
  dops2plot = find(plot_selection(1,5:9)==1);
  if any(dops2plot) & any(gps_dops),
    %%% PLOTPASS
    % generate the data for DOPS versus time
    n_plots = n_plots + 1;
    nrows = size(gps_dops,1);
    for k=1:length(dops2plot),
      dops_t((k-1)*nrows+1:k*nrows,1:2) = gps_dops(:,1:2);
      plotdops((k-1)*nrows+1:k*nrows,:) = [gps_dops(:,2+dops2plot(k))];
      prn_vis((k-1)*nrows+1:k*nrows,1:2) = [zeros(nrows,1) ones(nrows,1)*k]; 
      legend_save(k,:) = DOP_label(dops2plot(k),:);
    end;
    y_label = 'DOPS';
    title_str_dop = sprintf('DOPs for %s',label_string);
    fig_handles(n_plots) = plotpass(dops_t, plotdops, prn_vis, ...
                       title_str_dop, y_label);
    legend(legend_save,0)   
    set(fig_handles(n_plots),'position',[x_min-x_step*3 y_min x_max/2 y_max/2]);
  end; % if any(dops2plot) & any(t_dop_gps),

  if plot_selection(1,1)==1,

    % Create a figure for the sky plot and shift it left of elevation plot
    n_plots = n_plots + 1;
    fig_handles(n_plots) = figure('position',[x_min-x_step*4 y_min x_max/2 y_max/2], ...
      'NumberTitle','off','Color','Black','Units','pixels');
    set(fig_handles(n_plots),'Name',sprintf('Sky Plot for %s',label_string));

    % create a title for the sky plot
    t_plot_gps = visible_data(:,3:4);
    t_linear = gpst2sec(t_plot_gps);
    start_utc = gps2utc(t_plot_gps(1,:));
    title_str_sky = sprintf('%.1f Hours of Coverage from %d/%d/%d %d:%02d:%02d for %s',...
               (t_linear(length(t_linear)) - t_linear(1)) / 3600,...
               start_utc(2),start_utc(3),start_utc(1)-1900,start_utc(4:6), ...
               label_string);
    % save azel_gps az el t_los_gps prn_gps mask title_str_sky title_str_az title_str_el

    % generate a polar sky plot
    plot_handle = plotsky(visible_data(:,1), visible_data(:,2), ...
                  [visible_data(:,5), pass_numbers], title_str_sky);
  end;

%%%%% END ALGORITHM CODE %%%%%

% end of MAKEPLOT
