% ex_vect.m
%
% Function to demonstrate body-fixed masking, efficient breaking up of 
% the vectorized functions to conserve memory, changing GPS almanacs,
% and computing DGPS coverage and navigation solutions for an turning aircraft.
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

clear       % clear all variables in workspace
%close all  % close all open windows

% Set simulation parameters
d2r = pi/180;
mask_base = 0;             % simple 0 deg elevation Earth Mask (radians)
mask_b = [0 190 170;       % 0 elevation mask of fuselage between 190-170 deg
       30 170 190]*d2r;   % 30 deg elevation mask of tail between 170-190 deg 

% Input some receiver modeling parameters
ac_rec_code_noise = 3;       % 3 meters of code noise, typical mid quality C/A code
ac_rec_carrier_noise = .1;   % 10 cm of carrier noise       
ac_pr_err_model = [1 1 1 1 1 1];   % model all PR errors

sim_start_time_utc = [1997 9 17 0 0 0];   % Start Date (yr, mon, day, hr, mn, sec)
sim_stop_time_utc = [1997 9 17 0 5 0];    % Stop Date (yr, mon, day, hr, mn, sec)
sim_start_time_gps = utc2gps(sim_start_time_utc);  % Sim start in GPS time
sim_stop_time_gps = utc2gps(sim_stop_time_utc);    % Sim stop in GPS time
sim_start_time_sec = gpst2sec(sim_start_time_gps); % Sim start in GPS seconds
sim_stop_time_sec = gpst2sec(sim_stop_time_gps);   % Sim stop in GPS seconds

time_step = 1;          % Time step (sec)
alm_file = 'gps916.alm';        % GPS almanac file to be used here

% Set up the parameters to chunck up the run into multiple sections
% to add the capability of changing almanac/ephemeris data sets during 
% the simulation.  This allows for setting specific satellite healthy 
% and un-healthy during the simulation.  This same technique can be used 
% break up a simulation into time chunck to save RAM for very large simulations.
% See the example function for the large constellation broken into 
% time chunck while maintaining the vecotrization.

% Use Matlab 5 structures for readibility in almanac manipulation.  Start
% by developing 4 separate almanacs, all with the same satellite information.
almanac_data(1).alm = readyuma(alm_file);
almanac_data(2).alm = almanac_data(1).alm;
almanac_data(3).alm = almanac_data(1).alm;
almanac_data(4).alm = almanac_data(1).alm;

% Define start times for each almanac.  The start time will be used
% to determine when this almanac is switched into the simulation.  Beacuse
% this is a short simulation (300 seconds), almanac start times will be
% from the beginning of the sim (ie 0 -> 300 seconds).
almanac_data(1).start_time = 0;      % Start this alm at the beginning
almanac_data(2).start_time = 60;     % Start this alm at 1 minute
almanac_data(3).start_time = 240;    % Start this alm at 4 minutes
almanac_data(4).start_time = 250;    % Use this alm for the rest

% Set some satellites unhealthy for the different almanacs.  The satellite
% ID is in column #1 of the almanac, the health bit is coulmn #2.  The 
% value of 0 is healthy, anything else is unhealthy. For this example,
% we will set satellites 17 and 30 unhealty.
I_17 = find(almanac_data(1).alm(:,1) == 17);
I_30 = find(almanac_data(1).alm(:,1) == 30);
almanac_data(2).alm(I_17,2) = 1;    % set prn 17 unhealthy in almanac 2
almanac_data(3).alm(I_17,2) = 1;    % set prn 17 unhealthy in almanac 3
almanac_data(4).alm(I_30,2) = 1;    % set prn 30 unhealthy in almanac 4

% Read in airplane data:  time, position, velocity, attitude
load airplane.dat;                 % aircraft data in an ascii file
ac_time = airplane(:,1);           % aircraft relative time: 0-300 sec
ac_pos_all = [airplane(:,2:4)];    % aircraft ECEF xyz position (m)
ac_vel_all = [airplane(:,5:7)];    % aircraft ECEF xyz velocity (m/s)
ac_att_all = [airplane(:,8:10)];   % aircraft attitude wrt NED (deg)
ac_att_all = ac_att_all*d2r;       % convert to radians

% Convert aircraft relative time to gps time, using start_time as time0
ac_time = ac_time + sim_start_time_sec;   % absolute aircraft time (secs)
ac_time_s = sec2gpst(ac_time);            % GPS time version of ac_time

% Set the base station location to be at the aircraft starting location,
% but offset by about 1 km
base_loc_ecef = ac_pos_all(1,:) + [1000 0 0];

% Find out how many chuncks to break the simulation into.
num_chuncks = size(almanac_data,2);

%%%%% BEGIN ALGORITHM CODE %%%%%

% Initialize variables for storing the output for all the chuncks
t_nav_gps_all = [];
t_nav_dgps_all = [];
x_nav_gps_all = [];
x_nav_dgps_all = [];

t_vis_all = [];
t_dops_all = [];
dops_all = [];
num_sats_all = [];

% Loop over each data chunck
for ijk = 1:num_chuncks          

  % Set the start and stop time for this chunck
  this_start_sec = sim_start_time_sec + almanac_data(ijk).start_time;
  if ijk < num_chuncks
    this_stop_sec = sim_start_time_sec + ...
                    almanac_data(ijk+1).start_time - time_step;
  else
    this_stop_sec = sim_stop_time_sec;
  end % if ijk < num_chunck
  
  start_gps = sec2gpst(this_start_sec);
  stop_gps = sec2gpst(this_stop_sec);
  
  % Load the GPS almanac for the given almanac week
  alm_2_use = almanac_data(ijk).alm;

  % Sort out the unhealthy satellites
  I_gps_good = find(alm_2_use(:,2) == 0);
  alm_2_use = alm_2_use(I_gps_good,:);

  % Convert the almanacs to ephemeris format
  [gps_ephem] = alm2geph(alm_2_use);

  % Compute satellite positions in ECEF frame for the given time range and interval
  [t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem,start_gps,stop_gps,time_step);
  
  % Compute pseudo-range measurements for the base station
  [t_pr_base,prn_base,pr_base,pr_orb_base,pr_base_err] = ...
               pseudo_r(start_gps,base_loc_ecef,t_gps,[prn_gps x_gps],...
               mask_base,ac_pr_err_model,0);

  % Compute differential corrections at the base station
  [t_dpr, dpr] = ...
          diffcorr(t_pr_base, [prn_base pr_base], pr_orb_base, base_loc_ecef);

  % Compute pseudo-ranges for the aircraft receiver.
  % Set the masking to the base_mask value of 0 to block out satellites that
  % are below the horizon and therefore masked by the Earth.  It is important
  % to model the base station and remote PR at the same time in the simulation
  % and with the same seed.  Otherwise, the SA errors, because of the algorithms
  % used to generate them, will not be common.  If the SA errors are not common,
  % the DGPS solutions will not have any meaning. 
  [t_pr_ac,prn_ac,pr_ac,pr_orb_ac,pr_ac_err] = ...
               pseudo_r(ac_time_s, ac_pos_all, t_gps,[prn_gps x_gps],...
               mask_base,ac_pr_err_model,0,ac_rec_code_noise, ac_rec_carrier_noise);

  % Compute LOS vectors from airplane to GPS satellites, in ECEF
  [t_los, los_ac, ac_num, prn_los, ac_pos, gps_pos, I_ac, I_gps] = ...
        los(ac_time_s, ac_pos_all, t_gps, [prn_gps x_gps]);

  % Build up an attitude matrix that uses the input matrix, ac_att, but
  % is sync'd to the LOS vectors.  This is necessary since the LOS matrix
  % is 4-12 times the size of the ac_att matrix.  The LOS matrix has several
  % vectors at the same time point, and we need the attitudes to be correlated
  % with these LOS vectors.  A large attitude matrix, ac_att_large, is built
  % from ac_att using I_ac from losorbit.

  ac_att = ac_att_all(I_ac,:);  % make the large attitude matrix 

  % Convert LOS in ECEF to NED
  [los_ned] = ecef2ned(los_ac, ecef2lla(ac_pos));

  % Compute az and el in Earth frame
  [az_e el_e] = ned2azel(los_ned);

  % Find indices to satellites visible above the Earth mask
  I_vis_e = mask_vis(az_e, el_e, mask_base);

  % Reset the arrays to contain only data visible from the Earth
  if any(I_vis_e),
    t_los = t_los(I_vis_e,:);                 
    prn_los = prn_los(I_vis_e);
    gps_pos = gps_pos(I_vis_e,:);
    az_e = az_e(I_vis_e);
    el_e = el_e(I_vis_e);
    ac_att = ac_att(I_vis_e,:);
    los_ned = los_ned(I_vis_e,:);
  else
    fprintf('No visibile satellites during Earth masking.\n');
    return
  end;

  % At this point, all of the data vectors contain only data that has passed the
  % Earth visibility mask.

  % Rotate the Earth-visible LOS's into the airplane body frame
  los_body = ned2body(los_ned,ac_att); % Earth vis. LOSs in b frame

  % Convert LOS's in body frame to body-referenced az and el's
  % using ned2azel here since it is the same as body2azel.  This assumes that the
  % azimuth is defined relative to the body x-axis and elevation is relative
  % the the x-y body plane.
  [az_b,el_b] = ned2azel(los_body);

  % Apply the body-fixed masking model to the body az and el's
  I_vis_b = mask_vis(az_b,el_b,mask_b); % I_vis_b is the index of vis. sats in b

  % Use this index of visible satellites in the body frame to build the LOS data
  % set that is visible in the body frame, with all the maskings applied.  
  % Express this LOS data in the NED frame, however, to keep the DOPS and
  % Visibility calculations in the local-level frame.
  if any(I_vis_b),
    t_pr_ac = t_pr_ac(I_vis_b,:);     % PR times for the aircraft data
    prn_ac = prn_ac(I_vis_b);         % PR PRN for the aircraft data
    pr_ac = pr_ac(I_vis_b,:);         % PR for the aircraft data
    pr_orb_ac = pr_orb_ac(I_vis_b,:); % PR orbits for the aircraft data
    t_los = t_los(I_vis_b,:);         % times associated with LOS's
    los_ned = los_ned(I_vis_b,:);     % LOS vectors, in NED, visible in body
    los_body = los_body(I_vis_b,:);   % LOS vectors, in body frame, visible in body
    prn_los = prn_los(I_vis_b);       % visible GPS PRN
    gps_pos = gps_pos(I_vis_b,:);     % visible GPS orbit positions
  else
    fprintf('No visibile satellites during aircraft body masking.\n');
    return
  end;             

  % Compute a navigation solution with these PR measurements, no DGPS
  [t_nav_gps,x_nav_gps,num_dgps] = ...
     lsnav(t_pr_ac,pr_ac(:,1),[prn_ac pr_orb_ac],[ac_pos_all(1,:)+50 0]); 
  
  % Apply the differential corrections
  [t_cpr, prn_cpr, cpr, cpr_index] = ...
      add_dpr(t_pr_ac,[prn_ac pr_ac],t_dpr,dpr);
  
  % Compute a differential GPS (DGPS) navigation solution
  [t_nav_dgps,x_nav_dgps,num_dgps] = ...
     lsnav(t_cpr,cpr(:,1),[prn_ac(cpr_index) pr_orb_ac(cpr_index,:)],...
           [ac_pos_all(1,:) 0]);

  % Finally, compute DOPS and number of satellites tracked
  [dops,t_dops] = ned2dops(los_ned,t_los,-90*d2r);  % no masking
  
  [az_eb,el_eb] = ned2azel(los_body);        % az and el's for the body masking
  [t_vis, num_sats] = num_vis(az_eb, el_eb, t_los, -90*d2r);

  % Save the data for this time/almanac chunck
  t_nav_gps_all = [t_nav_gps_all; t_nav_gps];
  t_nav_dgps_all = [t_nav_dgps_all; t_nav_dgps];
  x_nav_gps_all = [x_nav_gps_all; x_nav_gps];
  x_nav_dgps_all = [x_nav_dgps_all; x_nav_dgps];
  
  t_vis_all = [t_vis_all; t_vis];
  num_sats_all = [num_sats_all; num_sats];
  t_dops_all = [t_dops_all; t_dops];
  dops_all = [dops_all; dops];
  
end % for ijk = 1:num_chuncks

%%% Plotting Section %%%

% generate the plot of visible satellites versus time
figure
plot_handle3 = plot_svn(t_vis_all, num_sats_all);
title('# of Visible Satellites');
ylabel('# of Visible Satellites');

% plot DOPS
figure
tm = gpst2sec(t_dops_all);       % time past GPS epoch (1980 in seconds)
tm = (tm - tm(1))/ 60;           % time past start in minutes
plot_handle = plot(tm,dops_all);
title('DOPS for Non-Differential Aircraft Navigation Solution')
xlabel('Time Past Aircraft Epoch (min)')
ylabel('DOPS')
legend('GDOP','PDOP','HDOP','VDOP','TDOP',0);                

% Plot GPS errors
gps_err = x_nav_gps_all(:,1:3) - ac_pos_all;
gps_err_ned = ecef2ned(gps_err,ecef2lla(ac_pos_all));
figure
tm = gpst2sec(t_nav_gps_all);
tm = (tm - tm(1))/60;
plot(tm, gps_err_ned);
title('GPS Errors for Aircraft Flight Profile')
xlabel('Time Past Aircraft Epoch (min)')
ylabel('meters')
legend('North','East','Down',0)

% Plot GPS errors
dgps_err = x_nav_dgps_all(:,1:3) - ac_pos_all;
dgps_err_ned = ecef2ned(dgps_err,ecef2lla(ac_pos_all));
figure
plot(tm, dgps_err_ned);
title('DGPS Errors for Aircraft Flight Profile')
xlabel('Time Past Aircraft Epoch (min)')
ylabel('meters')
legend('North','East','Down',0)

%%% End of plotting

% end of vis_ac.m
