function fig_handle = rdrxdata(auto_mode)

% Example function to generate a text box for displaying information about
% reading data from GPS receivers.
% Input:
%   auto_mode = optional flag (only include when in auto mode)

% Written by: Jimmy LaMance 9/6/97
% Copyright (c) 1998 by Constell, Inc.

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
   'Name','Reading GPS Receiver Data', ...
   'Tag','fign');
end;

fig_title_cell={'Reading GPS Receiver Data';};

descriptive_text_cell = ...
 {'For many applications, the simulation of GPS data is not sufficient and real';
  'data must be collected.  A common output format is the Nation Marine Electronics';
  'Association (NMEA) defined electrical and data protocol.  The data format is';
  'an ASCII based messaging system.  However, because of the variable message';
  'length and the fact that some messages span multiple records, this format is';
  'sometimes difficult to load into computer software.  The Constellation Toolbox provides';
  'a set of functions designed for reading and parsing NMEA messages.  Once the';
  'data is read into the Matlab environment, it can be saved in more convenient';
  'formats, such as ASCII files that are space or tab delimited, or the data can';
  'be analyzed using a combination of the Constellation Toolbox and built-in Matlab';
  'functions.';
  '';
  'Data that is available from the NMEA messages include position and velocity,';
  'Dilution of Precision (DOP) estimates from the receiver, satellites currently';
  'being used in the position solution with the signal to noise ratio, and the'; 
  'azimuth and elevation to all the visible satellites.';
  '';
  'For more information on the NMEA and the NMEA format, see the NMEA World';
  'Wide Web site at http://www4.coastalnet.com/nmea/default.html.';
  '   '};

if nargin < 1,
  textwin3(fig_title_cell,descriptive_text_cell,1);

  % Determine the location for the plot in upper right corner
  x_min = x_max / 10;
  y_min = y_max / 10;

  % set the figure position
  set(fig_handle,'position',[x_min y_min x_max/1.25 y_max/1.25]);
else,
  fig_handle = gcf;
  textwin3(fig_title_cell,descriptive_text_cell);
end;

% end of PLTDOPER
