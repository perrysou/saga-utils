function az_el_table = yuma2az_el(alm_file, location, start_time, stop_time, time_step)
% A function to read a YUMA file and output satellites' azimuth and
% elevation angles between start_time and stop_time for a given
% station lacation.
%%% This function is derived from iono_as_geomGen.m written by Jiyun Lee.
%  6/20/06 Godwin Zhang 
% 6/29/06 Seebany Datta-Barua Modified to keep information about 1024 week
% rollovers.  This file works with supertruth data for 10/30/03 and
% 301.alm,  This almanac has week rollover so that 10/30/03 is listed as week 218.
% 7 Aug 2013 SDB Commenting out the week rollover parts because the version
% of utc2gps.m I now use doesn't return that info. Output is now 5 cols:
% az_el_table = [mod(gps_week,1024) gps_sec prn az el].
%
% Output:
% table az_el_table = [mod(gps_week,1024) gps_sec rollover_flag prn az el]
% Inputs:   
%alm_file: YUMA filename string
%location:  The lla coordinate of a station, in radian.
%start_time and stop_time: Start and stop time in UTC format.
%time_step: time step in second
% This script uses Constellation Toolbox


% convert the station location from lat, long, alt. to ECEF vector
location_ecef = lla2ecef(location);

% load the GPS almanac for the given almanac week
alm_2_use = readyuma(alm_file);

%%% Don't need
% % % from Curt's comment : LAN = RAAN - GHA(Greenwich Hour Angle) 
% % % alm_2_use(:,8)  = alm_2_use(:,8) - 275.11*pi/180;%272.64*pi/180; (June 30,23:34:24)%275.11*pi/180; (June 27)% 279.06*pi/180; (July1) 

% sort out the unhealthy satellites
I_gps_good = find(alm_2_use(:,2) == 0);
alm_2_use = alm_2_use(I_gps_good,:);
I_gps_bad = find(alm_2_use(:,2) == ~0);
if ~isempty(I_gps_bad)
    disp('There are unhealthy satellites');
    I_gps_bad
end;


% convert the almanacs to ephemeris format
[gps_ephem] = alm2geph(alm_2_use);

% first convert the start and stop times to GPS time.
start_gps = utc2gps(start_time);
stop_gps = utc2gps(stop_time);

% compute satellite positions in ECEF frame for the given time range and interval
[t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem, start_gps, stop_gps, time_step);

% compute LOS vectors in ECEF frame
[t_los_gps, gps_los, los_ind] = los(t_gps(1,:), location_ecef, t_gps, [prn_gps x_gps]);

% convert LOS in ECEF to NED frame
[gps_los_ned] = ecef2ned(gps_los, location);

% Compute azimuth and elevation
[az, el] = ned2azel(gps_los_ned);

% Note: az_el_table will be nx6 because now t_gps is nx3, where the 3rd
% column is the rollover flag.
az_el_table = [t_gps, prn_gps, az, el];
% keyboard
% % first convert the start and stop times to GPS time.  Keep rollover info.
% [start_gps_week, start_gps_sec, start_day, start_rlvr_flag] = utc2gps(start_time);
% [stop_gps_week, stop_gps_sec, stop_day, stop_rlvr_flag] = utc2gps(stop_time);
% 
% % start_gps is 1x3, stop_gps is 1x3.  3rd column is the rollover flag.
% start_gps = [start_gps_week start_gps_sec start_rlvr_flag];
% stop_gps = [stop_gps_week stop_gps_sec stop_rlvr_flag];
% 
% % compute satellite positions in ECEF frame for the given time range and interval
% [t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem, start_gps, stop_gps, time_step);
% 
% % compute LOS vectors in ECEF frame
% [t_los_gps, gps_los, los_ind] = los(t_gps(1,:), location_ecef, t_gps, [prn_gps x_gps]);
% 
% % convert LOS in ECEF to NED frame
% [gps_los_ned] = ecef2ned(gps_los, location);
% 
% % Compute azimuth and elevation
% [az, el] = ned2azel(gps_los_ned);
% 
% % Note: az_el_table will be nx6 because now t_gps is nx3, where the 3rd
% % column is the rollover flag.
% az_el_table = [t_gps, prn_gps, az, el];
% % keyboard
