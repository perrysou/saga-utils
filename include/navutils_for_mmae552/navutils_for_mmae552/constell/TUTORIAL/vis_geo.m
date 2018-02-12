% vis_geo.m
%
% Function to demonstrate the visibility routines for GEO satellites. 
% Computes the azimuth and elevation to the GPS satellites and 
% displays an elevation vs. time, azimuth vs. time, # of satellites
% visible, and a sky plot showing azimuth/elevation pairs.
%
% See also vis_o, vis_e.m

% Written by: Jimmy LaMance
% Copyright (c) 1998 Constell, Inc.

% functions called: READYUMA, ALM2GEPH, UTC2GPS, PROPGEPH, LOS, ECEF2LL,
%                   NED2AZEL, VIS_DATA, PASSDATA, MAKEPLOT, KEP2GEPH, LL2DOPS

clear
close all

% Set simulation parameters

% Kepler Elements for the user satellite (not the GPS satellites)
% [sv_num, a (meters), e, i, node, argp, M, epoch week, epoch sec]
kep_elems = [1 42164172, 0, 0, 255*pi/180, 0*pi/180, 0*pi/180, 964, 0];

% Set the elevation mask for the GEO satellite.  In the satellite local
% level system, elevation is defined relative to the horizontal orbit plane
% with positive elvations in the "up" direction.  Therefore, this mask
% must be defined to be from -pi/2 up to 0 degrees elevation.  See the 
% description of the satellite local level system for more details (ECEF2LL.m)
mask_geo = [-pi/2 0 0 2*pi];            % [min_el max_el min_az max_az]

% Set the elevation mask for the GPS antenna (transmit).  The same local
% level coordinate system will apply. For an antenna on "top" of the GPS
% satellite a mask like the following would be used to describe an positive
% cone of elevation.
% mask_gps = 14.3*pi/180;
% See help on VIS_DATA for more details about entering masks and ECEF2LL
% for a description of the satellite local level system used to compute 
% the satellite masking.
gps_ant_cone = 45.3 * pi /180;                 % GPS antenna cone angle in rad
mask_gps = [-pi/2 -pi/2+gps_ant_cone 0 2*pi];  % [min_el max_el min_az max_az]

% Set the simulation start and stop time
start_time = [1998 7 4 0 0 0];     % Start Date (yr, mon, day, hr, mn, sc)
stop_time = [1998 7 4 1 0 0];     % Stop Date (yr, mon, day, hr, mn, sec)
time_step = 60;         % Time step (sec)
alm_file = 'gps916.alm';        % GPS almanac file

%%%%% BEGIN ALGORITHM CODE %%%%%

% load the GPS or GLONASS almanac for the given almanac week
alm_2_use = readyuma(alm_file);

% sort out the unhealthy satellites
I_gps_good = find(alm_2_use(:,2) == 0);
alm_2_use = alm_2_use(I_gps_good,:);

% convert the almanacs to ephemeris format
[gps_ephem] = alm2geph(alm_2_use);

% first convert the start and stop times to GPS time
start_gps = utc2gps(start_time);
stop_gps = utc2gps(stop_time);

% convert from the keplerian set to GPS ephemeris format
user_ephem = kep2geph(kep_elems);

% Compute user satellite positions in ECEF 
[t_user,prn_user,x_user,v_user] = ...
         propgeph(user_ephem,start_gps,stop_gps,time_step);

% Compute navigation satellite positions in ECEF 
[t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem,start_gps,stop_gps,time_step);

% Compute LOS vectors in ECEF 
[t_los_gps,los,los_ind,obscure_info] = los(t_user, x_user, t_gps, [prn_gps x_gps]);

% Rename some variables for ease of use later
gps_prn_los = prn_gps(los_ind(:,2));
I_user = los_ind(:,1);
I_gps = los_ind(:,2);

% Convert LOS in ECEF to local-level
[los_ll_geo] = ecef2ll(los, x_user(I_user,:),v_user(I_user,:));

% Convert LOS in ECEF to local-level (GPS) from the GPS to the GEO satellite.
% The LOS vector from the GPS satellite is the negative of the LOS vector
% to the GPS satellite.
[los_ll_gps] = ecef2ll(-los, x_gps(I_gps,:),v_gps(I_gps,:));

% Compute az and els 
[az_geo el_geo] = ned2azel(los_ll_geo);

% Find indices to GPS satellites above the GEO mask and not 
% obscured by the Earth
[az_el_geo, I_vis_geo] = vis_data(mask_geo, [az_geo, el_geo], obscure_info);
if any(I_vis_geo),
  % reset the arrays to contain only visible data
  t_los_gps = t_los_gps(I_vis_geo,:);
  los = los(I_vis_geo,:);
  los_ll_geo = los_ll_geo(I_vis_geo,:);
  los_ll_gps = los_ll_gps(I_vis_geo,:);
  az_geo = az_geo(I_vis_geo);
  el_geo = el_geo(I_vis_geo);
  gps_prn_los = gps_prn_los(I_vis_geo);
end;

% Compute az and els 
[az_gps el_gps] = ned2azel(los_ll_gps);

% Find indices to GPS satellites above the GEO mask and not 
% obscured by the Earth
[az_el_gps, I_vis_gps] = vis_data(mask_gps, [az_gps, el_gps]);
if any(I_vis_gps),
  % reset the arrays to contain only visible data
  t_los_gps = t_los_gps(I_vis_gps,:);
  los = los(I_vis_gps,:);
  los_ll_geo = los_ll_geo(I_vis_gps,:);
  az_geo = az_geo(I_vis_gps);
  el_geo = el_geo(I_vis_gps);
  gps_prn_los = gps_prn_los(I_vis_gps);
end;

% Compute DOPs
[dops, t_dops, num_sats] = ll2dops(los_ll_geo, t_los_gps);

%%% Plotting Section %%%
% Make arrays for using the MAKEPLOT function
vis_data = [az_geo el_geo t_los_gps gps_prn_los];
[pass_numbers, pass_times, pass_summary] = passdata(t_los_gps, 600, ...
         [ones(size(gps_prn_los)) gps_prn_los], vis_data(:,1:2));
num_vis = [t_dops,num_sats];
gps_dops = [t_dops dops];

% Compute pass information and print to screen
pass_times_utc = gps2utc(pass_times(:,2:3));
output_array = [pass_times_utc(:,2:3) pass_times_utc(:,1) pass_times_utc(:,4:6) ...
   pass_times(:,4)/60 pass_times(:,5:6) pass_summary(:,1,1)*180/pi ...
   pass_summary(:,2,1)*180/pi pass_summary(:,1,2)*180/pi ...
   pass_summary(:,2,2)*180/pi pass_summary(:,2,4)*180/pi];
fprintf('Start Time of Pass    Duration  Obs.  S/C Rise Az. Elev.  Set Az. Elev. Max Elev.\n');
fprintf('       (UTC)            (min)   PRN   PRN  (deg)   (deg)   (deg)  (deg)   (deg)\n');
fprintf('%2d/%2d/%4d %2d:%2d:%4.1f %7.2f %5d %5d %7.2f %6.2f %7.2f %6.2f %7.2f\n', ...
  output_array');

plot_selection = [0 1 1 1 1 1 1 1 1];
fig_handles = makeplot(vis_data, pass_numbers, num_vis, gps_dops,...
                       'GPS Visibility for a GEO Satellite.',plot_selection);

%%% End of plotting

% end of vis_o.m
