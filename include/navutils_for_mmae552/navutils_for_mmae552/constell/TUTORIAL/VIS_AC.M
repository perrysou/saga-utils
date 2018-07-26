% vis_ac.m
%
% Function to demonstrate body-fixed masking and Earth obscuring functions
% for an aircraft trajectory.
%
% An aircraft is modeled with an initial location at ECEF = (6378e3,0,0).  It
% has a ground speed of 200 m/s (450 mph).  The aircraft flies straight and
% level North for 60 sec, then banks to the right by 45 deg and executes a
% 360 deg. turn at constant altitude (about 128 sec).  At the end of the turn,
% the airplane levels back out and continues to fly North, straight and level, 
% until 300 secs have elapsed since the beginning of the flight.  Aircraft
% position, velocity, and attitude data are stored in a file at 1 sec intervals.
% This data is read in and used in this example.
%
% A GPS antenna is mounted to the top of the airplane.  It is masked by the
% fuselage (elevation < 0) and by the vertical tail (elevation < 30 deg for
% an azimuth range of 170 to 190 deg, relative to the nose of the airplane).
% A 5 deg mask of the Earth is also used.

% Written by: Jimmy LaMance
% Copyright (c) 1998 Constell, Inc.

% functions called: READYUMA, ALM2GEPH, UTC2GPS, GPST2SEC, SEC2GPST,
%                   PROPGEPH, LOS, ECEF2NED, NED2BODY, NED2AZEL, VIS_DATA,
%                   NED2DOPS, MAKEPLOT

clear       % clear all variables in workspace
close all   % close all open windows

% Set simulation parameters
d2r = pi/180;

% Define the body fixed antenna mask
mask_b = [0 190 170;      % 0 elevation mask of fuselage between 190-170 deg
         30 170 190]*d2r; % 30 deg elevation mask of tail between 170-190 deg

% Define an absolute start time for the simulation.  The aircraft trajectory
% is given in relative time (seconds past the start of the trajectory).  Get
% the start time in seconds past the GPS epoch.
start_gps = [916 0];
start_sec = gpst2sec(start_gps);

% Define the GPS almanac file to use
alm_file = 'gps916.alm';         % GPS almanac file to be used here

%%%%% BEGIN ALGORITHM CODE %%%%%

% Load the GPS almanac for the given almanac week
alm_2_use = readyuma(alm_file,0);

% Sort out the unhealthy satellites
I_gps_good = find(alm_2_use(:,2) == 0);
alm_2_use = alm_2_use(I_gps_good,:);

% Convert the almanacs to ephemeris format
[gps_ephem] = alm2geph(alm_2_use);

% Read in airplane data:  time, position, velocity, attitude
load airplane.dat;              % aircraft data in an ascii file
ac_time = airplane(:,1);        % aircraft relative time: 0-300 sec
ac_pos = [airplane(:,2:4)];     % aircraft ECEF xyz position (m)
ac_vel = [airplane(:,5:7)];     % aircraft ECEF xyz velocity (m/s)
ac_att = [airplane(:,8:10)];    % aircraft attitude wrt NED (deg)
ac_att = ac_att*d2r;            % convert to radians

% Convert aircraft relative time to gps time, using start_time as time0
ac_time = ac_time + start_sec;  % absolute aircraft time (secs)
ac_time_s = sec2gpst(ac_time);  % GPS time version of ac_time

% Compute satellite positions in ECEF frame at the times from the aircraft
% trajectory.  This is a simple way to obtain time synced GPS positions and
% velocities if the desired output times are not in even steps.
[t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem,ac_time_s);

% Compute LOS vectors from airplane to GPS satellites, in ECEF. The los_ind 
% (line-of-site indices) will also be used to line up all of the data with 
% the positions used to compute the LOS vector.   In addition, we will get the 
% Earth obscuring information back from LOS to determine which satellites
% are blocked by the Earth.  Get the tangent altitude 
% (how far above the surface of the Earth does the LOS vector pass) relative
% to the spheriod.  
earth_model = 0;           % 1 = use ellipsoid, 0 = use spheriod
[t_los_gps, los_ecef, los_ind, obscure_info] = ...
     los(ac_time_s, ac_pos, t_gps, [prn_gps x_gps],earth_model);

% Rename some of the variables for easier use later.
prn_nav = prn_gps(los_ind(:,2));
gps_pos = x_gps(los_ind(:,2),:);
I_ac = los_ind(:,1);
  
% Build up an attitude matrix that uses the input matrix, ac_att, but
% is sync'd to the LOS vectors.  This is necessary since the LOS matrix
% is 4-12 times the size of the ac_att matrix.  The LOS matrix has several
% vectors at the same time point, and we need the attitudes to be correlated
% with these LOS vectors.  A large attitude matrix, ac_att_large, is built
% from ac_att using I_ac from losorbit.
ac_att_large = ac_att(I_ac,:);      % make the large attitude matrix 
ac_pos_large = ac_pos(I_ac,:);  % do the same for the airplane position data

% Rotate the LOS's from ECEF to the NED frame
los_ned = ecef2ned(los_ecef,ecef2lla(ac_pos_large));

% Rotate the Earth-visible LOS's into the airplane body frame
los_body = ned2body(los_ned,ac_att_large); 

% Convert LOS's in body frame to body-referenced az and el's
% (use ned2azel here since it is the same as body2azel would be, if we had one)
[az_b,el_b] = ned2azel(los_body);

% Apply the body-fixed masking model to the body az and el's.
% Use the obscure_info from LOS to mask out the satellite that are blocked by
% the Earth.  To demonstrate the use of this capability, we have chosen to
% mask all satellites that the LOS would be within 5 km of the Earth using
% a spherical Earth model.  This can also be done with an elliptical 
% Earth model.
obscure_height = 5000;    % min height above the spheriod/ellipse (m) 
[az_el, I_vis_b] = vis_data(mask_b,[az_b,el_b],obscure_info,obscure_height);    

% Use this index of visible satellites in the body frame to build the LOS data
% set that is visible in the body frame, with all the maskings applied.  
% Express this LOS data in the NED frame, however, to keep the DOPS and
% visibility calculations in the local-level frame.
if any(I_vis_b),
  los_ned = los_ned(I_vis_b,:);      % LOS vectors, visible, in NED frame
  los_body = los_body(I_vis_b,:);    % LOS vectors, visible, in body frame
  t_los_gps = t_los_gps(I_vis_b,:);  % times associated with above LOS's
  prn_nav = prn_nav(I_vis_b);        % prn's to match the LOS's
end;    

% Finally, compute DOPS, number of satellites tracked, and the NED azimuth
% and elevation for plotting.
[dops,t_dops,num_sats_vis] = ned2dops(los_ned,t_los_gps);
[az_ned,el_ned] = ned2azel(los_ned);

% Make arrays for using the MAKEPLOT function
ac_vis_data = [az_ned el_ned t_los_gps prn_nav];
pass_numbers = ones(size(ac_vis_data,1),1);
num_ac_vis = [t_dops num_sats_vis];
gps_dops = [t_dops dops];

fig_handles = makeplot(ac_vis_data, pass_numbers, num_ac_vis, gps_dops,...
                       'Aircraft Example for GPS Visibility with Attitude');

% end of vis_ac.m
