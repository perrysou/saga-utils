% exdgpsac.m
%
% Function to demonstrate body-fixed masking, efficient breaking up of 
% the vectorized functions to conserve memory, changing GPS almanacs,
% and computing DGPS coverage and navigation solutions for an aircraft
% which models aircraft atitude and it's effects on DGPS satellite visibility. 
%
% This example has a large body-fixed masking (20 degress over about 1/2 of the
% sky) to demonstrate the effects of bad satellite geometry.  The effects
% on the differentially corrected solutions and the non-differentially 
% corrected solution is apparent.  DGPS errors in excess of 20 meters
% and raw GPS errors in excess of 600 meters are shown.  These large
% errors result from bad geometry.  To test this, set the aircraft masking 
% to -10 deg over the entire azimuth and see the improvement. However, because
% of the aircraft maneuvers and the body-fixed masking, the satellite geometry
% (observed as DOP values) is still poor.
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
close all  % close all open windows

% Set simulation parameters
d2r = pi/180;
mask_base = 0;             % simple 0 deg elevation Earth Mask (radians)
mask_b = [0 190 170;       % 0 elevation mask of fuselage between 190-170 deg
         20 170 190]*d2r;  % 30 deg elevation mask of tail between 170-190 deg 

% Input some receiver modeling parameters
ac_rec_code_noise = 3;       % 3 meters of code noise, typical mid quality C/A code
ac_rec_carrier_noise = .1;   % 10 cm of carrier noise       

% Set up the PR modeling for the aircraft.  Model SA epsilon, SA dither,
% troposphere, ionosphere, receiver clock, and receiver white noise (1's).
% Do not model line bias, Earth rotation, satellite motion, 
% satellite clocks, or relativity.  This set of modeling will be used to 
% simulate PR measurements for both the base station and the aircraft.
ac_pr_err_model = [1 1 1 1 1 1 0 0 0 0 0];   % model all PR errors

% Set the seed value to use in generate the PR model errors
seed = 0;

sim_start_time_utc = [1997 9 17 0 0 0];   % Start Date (yr, mon, day, hr, mn, sec)
sim_stop_time_utc = [1997 9 17 0 5 0];    % Stop Date (yr, mon, day, hr, mn, sec)
sim_start_time_gps = utc2gps(sim_start_time_utc);  % Sim start in GPS time
sim_stop_time_gps = utc2gps(sim_stop_time_utc);    % Sim stop in GPS time
sim_start_time_sec = gpst2sec(sim_start_time_gps); % Sim start in GPS seconds
sim_stop_time_sec = gpst2sec(sim_stop_time_gps);   % Sim stop in GPS seconds

time_step = 1;          % Time step (sec)
alm_file = 'gps916.alm';        % GPS almanac file to be used here

% Set up the parameters to chunk up the run into multiple sections
% to add the capability of changing almanac/ephemeris data sets during 
% the simulation.  This allows for setting specific satellite healthy 
% and un-healthy during the simulation.  This same technique can be used 
% break up a simulation into time chunk to save RAM for very large simulations.
% See the example function for the large constellation broken into 
% time chunk while maintaining the vecotrization.

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
% but offset by 1 km in the X-ECEF direction.
base_loc_ecef = ac_pos_all(1,:) + [1000 0 0];  
base_vel = zeros(size(base_loc_ecef));

% Find out how many chunks to break the simulation into.
num_chunks = size(almanac_data,2);

%%%%% BEGIN ALGORITHM CODE %%%%%

% Initialize variables for storing the output for all the chunks
t_nav_gps_all = [];
t_nav_dgps_all = [];
x_nav_gps_all = [];
x_nav_dgps_all = [];

t_vis_all = [];
t_dops_all = [];
dops_all = [];
num_sats_all = [];

% Loop over each data chunk
for ijk = 1:num_chunks          

  % Set the start and stop time for this chunk
  this_start_sec = sim_start_time_sec + almanac_data(ijk).start_time;
  if ijk < num_chunks
    this_stop_sec = sim_start_time_sec + ...
                    almanac_data(ijk+1).start_time - time_step;
  else
    this_stop_sec = sim_stop_time_sec;
  end % if ijk < num_chunk
  
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
  
  % Compute pseudo-range measurements for the base station.  Get the indicies
  % used to compute the PR (which base station position goes with which 
  % GPS satellite position) and the Earth obscure information.  
  [t_pr_base,prn_base,pr_base,pr_base_err,prn_ndex] = ...
               pseudo_r(start_gps,base_loc_ecef, base_vel, ...
               t_gps,[prn_gps x_gps],v_gps,ac_pr_err_model,seed);

  % Use the indices from PSEUDO_R to get the base station and GPS satellite
  % positions for use in creating the LOS.  This will be used to compute
  % satellite masking.
  base_los = x_gps(prn_ndex(:,2),:) - base_loc_ecef(prn_ndex(:,1),:);
  
  % Convert the ECEF LOS to NED
  base_los_ned = ecef2ned(base_los,ecef2lla(base_loc_ecef,1));

  % Compute azimuth and elevation from the NED vectors
  [az, el] = ned2azel(base_los_ned);
    
  % Compute masking in the NED frame.  Since this is a base station, no obscure
  % information needs to be used.  Instead a simple az/el masking for the
  % antenna pattern is modeled.
  [visible_data, I_pass] = vis_data(mask_base, [az el]);
  
  % Eliminate the PR data not passing the masking test
  t_pr_base = t_pr_base(I_pass,:);
  prn_base = prn_base(I_pass,:);
  pr_base = pr_base(I_pass,:);
  pr_base_err = pr_base_err(I_pass,:);
  prn_ndex = prn_ndex(I_pass,:);
  
  % Compute differential corrections at the base station.  Do not use any
  % pseudo-range models (SA, atmospheric, satellite motion, etc.).  All 
  % of these parameters that are modeled in the raw PR measurements will be
  % effectively removed from the base station PR with the differential 
  % corrections.
  [t_dpr, dpr] = diffcorr(t_pr_base, [prn_base pr_base], ...
                          x_gps(prn_ndex(:,2),:), v_gps(prn_ndex(:,2),:),...
                          base_loc_ecef);

  % Compute pseudo-ranges for the aircraft receiver.  Use increased code and
  % carrier noise to more closely imitate an aircraft type receiver.
  % Model the base station and remote PR at the same time in the simulation
  % and with the same seed.  Otherwise, the SA errors, because of the algorithms
  % used to generate them, will not be common.  If the SA errors are not common,
  % the DGPS solutions will not have any meaning. Also get the Earth obscuring
  % information back from PSEUDO_R.  This will allow the Earth to be masked
  % out with VIS_DATA.  The default Earth obscuring model is a spherical Earth.
  [t_pr_ac,prn_ac,pr_ac,pr_ac_err,prn_ndex,obscure_info] = ...
               pseudo_r(ac_time_s, ac_pos_all, ac_vel_all, ...
               t_gps,[prn_gps x_gps],v_gps,...
               ac_pr_err_model,seed,ac_rec_code_noise, ...
               ac_rec_carrier_noise);

  % Rename some variables for ease of use later
  prn_los = prn_gps(prn_ndex(:,2));
  ac_pos = ac_pos_all(prn_ndex(:,1),:);
  gps_pos = x_gps(prn_ndex(:,2),:);
  I_ac = prn_ndex(:,1);
  I_gps = prn_ndex(:,2);
  los_ac = gps_pos - ac_pos;
  
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

  % Rotate the Earth-visible LOS's into the airplane body frame
  los_body = ned2body(los_ned,ac_att); % Earth vis. LOSs in b frame

  % Convert LOS's in body frame to body-referenced az and el's
  % using ned2azel here since it is the same as body2azel.  This assumes that the
  % azimuth is defined relative to the body x-axis and elevation is relative
  % the the x-y body plane.
  [az_b,el_b] = ned2azel(los_body);

  % Apply the body-fixed masking model to the body az and el's.  Use the 
  % obscure information from the PSEUDO_R routine to eliminate data that is
  % masked by the Earth.  The default (no input for minimum tangent altitude)
  % is at the surface of the Earth.  See vis_ac for using a minimum tangent
  % altitude with the obscure ifnromation.
  [azel_vis, I_vis_b] = vis_data(mask_b,[az_b el_b],obscure_info); 

  % Use this index of visible satellites in the body frame to build the LOS data
  % set that is visible in the body frame, with all the maskings applied.  
  % Express this LOS data in the NED frame, however, to keep the DOPS and
  % Visibility calculations in the local-level frame.
  if any(I_vis_b),
    t_pr_ac = t_pr_ac(I_vis_b,:);     % PR times for the aircraft data
    prn_ac = prn_ac(I_vis_b);         % PR PRN for the aircraft data
    pr_ac = pr_ac(I_vis_b,:);         % PR for the aircraft data
    pr_orb_ac = gps_pos(I_vis_b,:);   % PR orbits for the aircraft data
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
  [dops,t_dops,num_sats] = ned2dops(los_ned,t_pr_ac);         % no masking
  
  [az_eb,el_eb] = ned2azel(los_body);        % az and el's for the body masking

  % Save the data for this time/almanac chunk
  t_nav_gps_all = [t_nav_gps_all; t_nav_gps];
  t_nav_dgps_all = [t_nav_dgps_all; t_nav_dgps];
  x_nav_gps_all = [x_nav_gps_all; x_nav_gps];
  x_nav_dgps_all = [x_nav_dgps_all; x_nav_dgps];
  
  % Store the data to be stored for all of the chunks.  These variables are 
  % initialized before the start of the chunking loop.  This is where you add
  % additional variable to be stored for output.  Be sure to initialize
  % the same way these are.
  t_vis_all = [t_vis_all; t_dops];
  num_sats_all = [num_sats_all; num_sats];
  t_dops_all = [t_dops_all; t_dops];
  dops_all = [dops_all; dops];
                           
  % Writing out data to a file...
  % To write all of the data out to a file uncomment the following code
     %write_string = sprintf('save exdata%02d.mat',ijk);
     %eval(write_string);
                             
  % To write a part of the data out to a file uncomment the following code
  % and add the variables you're interested in.
      %write_string = sprintf('save exdata%02d.mat t_vis_all dops_all',ijk);
      %eval(write_string);            
  
  % To read in this data and put it back together, use the load command
  % and the same concatenation procedure as above to generate matrices with
  % all of the data.

end % for ijk = 1:num_chunks

%%% Plotting Section %%%
% generate the plot of visible satellites versus time
fh(1) = plotpass(t_vis_all,num_sats_all,ones(size(num_sats_all)),...
                    'Number of Common View (Base & A/C) Satellites', '# Visible');

% plot DOPS
fh(2) = figure;
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
fh(3) = figure;
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
fh(4) = figure;
plot(tm, dgps_err_ned);
title('DGPS Errors for Aircraft Flight Profile')
xlabel('Time Past Aircraft Epoch (min)')
ylabel('meters')
legend('North','East','Down',0)

% Stack the figures such that they are not on top of each other
x_wide = .6;
y_high = .6;
x_start = .03;
y_loc = .08;
x_step = .1;

for i = 1:length(fh)
  x_loc = x_start + (i-1) * x_step;
  pos = [x_loc y_loc x_wide y_high];
  set(fh(i),'Units','normalized');
  set(fh(i),'Position',pos);
end % for i = 1:legnth(fh)

%%% End of plotting

% end of exdgpsac.m
