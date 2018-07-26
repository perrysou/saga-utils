function fig_handle = rdyumtxt(auto_mode)

% Example function to diplay text about YUMA formatted almanacs.
% Input:
%   auto_mode = optional flag (only include when in auto mode)

% Written by: Jimmy LaMance 9/5/97
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
   'Name','Reading YUMA Almanac Files', ...
   'Tag','fign');
end;

fig_title_cell={'Reading YUMA Almanac Files';};

descriptive_text_cell = ...
 {'The almanac data for the GPS and GLONASS constellations provides';
  'information on the orbits for those satellites.  For GPS satellites, this';
  'information is broadcast in the subframe 4 and 5 messages.  Each GPS satellite';
  'broadcasts the almanac data for all of the satellites.  This allows a receiver';
  'to predict the locations for the entire constellation of satellites by decoding';
  'the data from a single satellite.  The Constellation Toolbox utilizes the almanac data';
  'to compute satellite visibility, generate satellite positions for simulations,';
  'and examine the current health of the constellation.  This gives a user access';
  'to current and past satellite constellation information.';
  '';
  'The YUMA format is an ASCII standard for the almanac data.  Almanac data';
  'in YUMA format is available from many commercial receivers.  Current as well as';
  'past YUMA formatted almanac data is available from the U.S. Coast Guard at their';
  'World Wide Web site http://www.navcen.uscg.mil.';
  ''; 
  'The almanac data is a subset of the more precise ephemeris data.  To compute';
  'satellite orbits based on the almanac data, the ephemeris information not provided';
  'in the almanac is filled with appropriate data based on the satellite orbit. This';
  'creates a pseudo-ephemeris which is used with the standard ICD-GPS 200';
  'ephemeris propagator to predict satellite positions.'};

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
