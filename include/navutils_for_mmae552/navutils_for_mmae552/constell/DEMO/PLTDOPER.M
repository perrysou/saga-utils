function fig_handle = pltdoper(auto_mode)

% Example function to plot processed DOP contribution to the position error.
% The data are loaded from the file demoplt1.mat
% which is generated by the example function demoplt1.m
% Input:
%   auto_mode = optional flag (only include when in auto mode)

% Written by: Jimmy LaMance 9/5/97
% Copyright (c) 1998 by Constell, Inc.

% load in the data
load demoplt1;

if nargin < 1,
  % Get the screen size in pixels to use for location of plots
  set(0,'units','pixels');
  screen_size = get(0,'screensize');
  y_max = screen_size(2) + screen_size(4) - 60;
  x_max = screen_size(1) + screen_size(3) - 50;
  x_step = 110;
  y_step = 60 + y_max;

  % set the figure colors to be black background like Matlab 4
  colordef none;                               

  % generate the position error figure
  fig_handle = figure('color','black', ...
   'position',[30 50 x_max-30 y_max-50],  ...
   'NumberTitle','off', ...
   'Name','DOP Contribution to Position Error', ...
   'Tag','fign');
end;

fig_title_cell={'DOP Contribution to Position Error';};
x = t_dop_plot ./ 60;
y = pos_err_dops;
axis_label_cell = {'Minutes Since Start'; 'Position Error (m)'; ...
  'Estimated Position Errors Computed from DOPs';};

legend_cell = {'PDOP Error';...
               'HDOP Error';...
               'VDOP Error'};

descriptive_text_cell = ...
 {'Dilution of Precision (DOP) is a measure of the ability of the satellite';
  'geometry to provide a position or velocity solution.  Lower values of DOP';
  'indicate a stronger geometry, leading to better position solutions.  A good';
  'first estimate of position (or velocity) accuracy is to multiply the DOP value by the';
  'uncertainty or sigma of the pseudorange (PR) or Doppler measurement.  Shown here ';
  'are the position error estimates based on the DOPs for a PR sigma of 23 meters.'};

if nargin < 1,
  text_win(fig_title_cell,x,y,axis_label_cell,legend_cell, ...
      descriptive_text_cell,1);

  % Determine the location for the plot in upper right corner
  x_min = x_max / 10;
  y_min = y_max / 10;

  % set the figure position
  set(fig_handle,'position',[x_min y_min x_max/1.25 y_max/1.25]);
else,
  fig_handle = gcf;
  text_win(fig_title_cell,x,y,axis_label_cell,legend_cell, ...
      descriptive_text_cell);
end;

% end of PLTDOPER
