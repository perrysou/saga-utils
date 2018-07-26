% vis_o.m
%
% Function to demonstrate the visibility routines for orbiting satellites. 
% Computes the azimuth and elevation to the GPS satellites and 
% displays an elevation vs. time, azimuth vs. time, # of satellites
% visible, and a sky plot showing azimuth/elevation pairs.
%
% See also vis_e.m

% Written by: Maria Evans
% Copyright (c) 1998 Constell, Inc.

% functions called: READYUMA, ALM2GEPH, UTC2GPS, PROPGEPH, LOS, ECEF2LL,
%                   NED2AZEL, VIS_DATA, PASSDATA, MAKEPLOT, KEP2GEPH, LL2DOPS

clear
close all

% Set simulation parameters

% Kepler Elements for the user satellite (not the GPS satellites)
% [sv_num, a (meters), e, i, node, argp, M, epoch week, epoch sec]
kep_elems = [1 7000000, 0, 45*pi/180, 0*pi/180, 0*pi/180, 0*pi/180, 900, 0];

mask = 0.0*pi/180;          % 0 radians elevation Mask
start_time = [1997 9 17 0 0 0];     % Start Date (yr, mon, day, hr, mn, sc)
stop_time = [1997 9 17 0 30 0];         % Stop Date (yr, mon, day, hr, mn, sec)
time_step = 60;         % Time step (sec)
alm_file = 'gps916.alm';        % GPS almanac file

%%%%% BEGIN ALGORITHM CODE %%%%%

% load the GPS or GLONASS almanac for the given almanac week
alm_2_use = readyuma(alm_file,0);

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

% convert LOS in ECEF to local-level
[los_ll] = ecef2ll(los, x_user(I_user,:),v_user(I_user,:));

% compute az and els 
[az el] = ned2azel(los_ll);

% Find indices to satellites above the mask and not obscured by the Earth
[az_el, I_vis] = vis_data(mask, [az, el], obscure_info);
if any(I_vis),
  % reset the arrays to contain only visible data
  t_los_gps = t_los_gps(I_vis,:);
  los_ll = los_ll(I_vis,:);
  az = az(I_vis);
  el = el(I_vis);
  gps_prn_los = gps_prn_los(I_vis);
end;

% Compute DOPs
[dops, t_dops, num_sats] = ll2dops(los_ll, t_los_gps);

%%% Plotting Section %%%
% Make arrays for using the MAKEPLOT function
vis_data = [az el t_los_gps gps_prn_los];
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

fig_handles = makeplot(vis_data, pass_numbers, num_vis, gps_dops,...
                       'GPS Visibility from a=7000 km, i=45 deg.');

%%% End of plotting

% end of vis_o.m
