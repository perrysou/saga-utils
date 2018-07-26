% ex_gps.m
%
% Example script for generating GPS measurements and computing position
% and velocity solutions using the raw data. This example is for a static 
% receiver fixed to the surface of the Earth.

% Written by: Jimmy LaMance 
% Copyright (c) 1998 by Constell, Inc.

% functions called: UTC2GPS, FIND_ALM, READYUMA, ALM2GEPH,
%                   PROPGEPH, ECEF2NED, PSEUDO_R, ECEF2LLA,
%                   NED2AZEL, VIS_DATA, DOPPLER, LSNAV,
%                   LOS, NED2DOPS, DOPS2ERR, PLOTPASS

% Variables to be used in the examples
start_time = [1998 7 4 0 0 0];         % UTC start time
duration = 1*3600;                     % 1 hour duration
base_station = [40*pi/180 255*pi/180 2000];       % base station coordinates

sample_rate = 20;                   % data sampling every 20 seconds

% Conversion from degree to radians
d2r = pi / 180;

% Compute the stop time based on the start time and duration
start_gps = utc2gps(start_time);
start_gps_sec = gpst2sec(start_gps);
stop_gps_sec = start_gps_sec + duration;
stop_gps = sec2gpst(stop_gps_sec);

% Find the almanac that is most recent to the start time
alm_file = find_alm(start_gps(1));

% Read in the almanac
alm = readyuma(alm_file,0);

% Sort out the unhealthy satellites
I_healthy = find(alm(:,2) == 0);
alm = alm(I_healthy,:);

% Convert from almanac to ephemeris format
gps_ephem = alm2geph(alm);

% Compute satellite positions in ECEF 
[t_gps,prn_gps,x_gps,v_gps] = propgeph(gps_ephem,start_gps,stop_gps,sample_rate);

% Set the remote time matrix to be the same as the base station
t_base = t_gps(1,:);

% Convert the base and remote station coordinates to ECEF
x_base = ones(size(t_base,1),1) * lla2ecef(base_station);  
v_base = zeros(size(x_base));

% Generate PR measurements for the base and remote stations using default
% values for masking, SA modeling, and random number seeding
% for the base station, force the modeling to have no receiver clock error
% to simulate the base station capability                                       
base_mask = 10*d2r;

% Model SA, tropo, iono, and receiver clock.  Do not model satellite motion,
% Earth rotation, satellite clocks, line baises, or relativity
base_model = [1 1 1 1 1 1 0 0 0 0 0];   % turn all all error models
base_code_noise = .2;
base_carrier_noise = .01;
base_seed = 0; 
max_latency = 2;            % set the maximum latency for a differential corr.

[t_pr_base,prn_base,pr_base,pr_errors_base,base_ndex] = ...
               pseudo_r(t_base,x_base,v_base,t_gps,[prn_gps x_gps],...
                        v_gps,base_model,base_seed, ...
                        base_code_noise, base_carrier_noise);

% Compute Doppler measurements at the site to be used 
% in computing velocity solutions
doppler_model = [1 1 1]; 
dop_noise = .3;
[t_dop,prn_dop,dopp,dop_orb,dop_err] = ...
     doppler(t_base,x_base,v_base,t_gps,[prn_gps x_gps],v_gps,...
             doppler_model,base_seed,dop_noise);

% Compute masking for the base station.  Start by compute the LOS from
% the base to the satellites in ECEF.  This is easily done using the indices
% returned from PSEUDO_R where the LOS has been vectorized.  Masking is
% done external to PSEUDO_R so that it can be performed in any coordinate
% system.  This example is for a station on the surface of the Earth with
% the antenna boresight pointed up.  Therefore, the NED system is used.
los_base_ecef = x_gps(base_ndex(:,2),:) - x_base(base_ndex(:,1),:);

% Rotate the LOS to NED, using the base station as the reference for the NED
% coordinate system and rotation.
ref_lla = ecef2lla(x_base(base_ndex(:,1),:));
los_base_ned = ecef2ned(los_base_ecef,ref_lla);

% Compute the azimuth/elevation of the NED vectors
[az, el] = ned2azel(los_base_ned);

% Apply the masking in the NED coordinate system
[visible_data, I_pass] = vis_data(base_mask,[az el]);

% Remove all of the base station data that did not pass the masking test
t_pr_base = t_pr_base(I_pass,:);
prn_base = prn_base(I_pass,:);
pr_base = pr_base(I_pass,:);   
dopp = dopp(I_pass,:);
dop_orb = dop_orb(I_pass,:);
pr_errors_base = pr_errors_base(I_pass,:);
base_ndex = base_ndex(I_pass,:);
x_ned = los_base_ned(I_pass,:);

% Rename the pr_errors for easier use later
pr_sa_errb = pr_errors_base(:,1) + pr_errors_base(:,2);
trop_base = pr_errors_base(:,3) + pr_errors_base(:,4);
iono_base = pr_errors_base(:,5);
clk_biasb = pr_errors_base(:,6);
clk_driftb = pr_errors_base(:,7);

% Compute a position and velocity solution at the remote site without DGPS
[t_nav,x_nav,num_remote,nav_index,v_nav] = ...
     lsnav(t_pr_base,pr_base(:,1),...
           [prn_base x_gps(base_ndex(:,2),:)],[x_base(1,:) 0],dopp, ...
            dop_orb(:,4:6),[0 0 0 0]); 

% Compute DOPs at the remote station (with masking for the remote station)
[dops, t_dops, num_base] = ned2dops(x_ned,t_pr_base);

% Estimate position errors from DOPs
sigma_pr = 40;
[pos_err_dops] = dops2err(dops,sigma_pr);
     
% Compute the position errors with and without DGPS
% these vectors are in ECEF
pos_err = x_nav(:,1:3) - ones(size(x_nav,1),1) * x_base;
vel_err = v_nav(:,1:3) - ones(size(v_nav,1),1) * v_base;

% Rotate position difference to NED relative to the truth location
% for easier interpretation                                       
[pos_err_ned] = ecef2ned(pos_err, ecef2lla(x_base));
[vel_err_ned] = ecef2ned(vel_err, ecef2lla(x_base));

% Generate a figure with 2 plots, one with the uncorrected position errors
% and one with differentially corrected position solutions     
plot_time = t_nav(:,2) - t_nav(1,2);  % simple time in seconds

% generate the position error figure
clear fh
fh(1) = figure;
subplot(2,1,1)
plot(plot_time,pos_err_ned)
legend('North','East','Down',-1);
ylabel('Position Err (m)')                
title_string = ...
  sprintf('Position Errors without DGPS for a Static Receiver');
title(title_string) 

subplot(2,1,2)
plot(vel_err_ned)
legend('North','East','Down',-1);
ylabel('Velocity Err (m/s)')
title_string = ...
  sprintf('Velocity Errors without DGPS for a Static Receiver');
title(title_string) 
xlabel('time past start (sec)') 

% generate the plots for S/A errors
fh(2) = plotpass(t_pr_base,pr_sa_errb,prn_base,...
         'Example of SA Errors on Pseudorange Measurements',...
         'SA PR err (m)');
         
% plot troposphere contribution to the PR error
fh(3) = plotpass(t_pr_base,trop_base,prn_base,...
         'Example of Troposphere Effects on Pseudorange Measurements',...
         'Tropo (m)');

% plot ionosphere contribution to the PR error
fh(4) = plotpass(t_pr_base,iono_base,prn_base,...
         'Example of Ionosphere Effects on Pseudorange Measurements',...
         'Iono (m)');

% plot the receiver clock bias for a single satellite 
% the receiver clock bias is common to all satellites
I = find(prn_base == prn_base(2));
fh(5) = plotpass(t_pr_base(I,:),clk_biasb(I),prn_base(I),...
         'Example of Receiver Clock Bias on Pseudorange Measurements',...
         'Clock Bias (m)');

fh(6) = plotpass(t_pr_base(I,:),clk_driftb(I),prn_base(I),...
         'Example of Receiver Clock Drift on Pseudorange Measurements',...
         'Clock Drift (m/2)');

% generate a plot of the resulting position error from the DOPs
t_dop_plot = (t_dops(:,1) * 604800 + t_dops(:,2)) - ...
             (t_dops(1,1) * 604800 + t_dops(1,2));

fh(7) = figure;
plot(t_dop_plot,pos_err_dops);
legend('PDOP Error','HDOP Error','VDOP Error',0);
ylabel('Position Err (m)')
title('Estimated Position Errors Computed from DOPs') 
xlabel('time past start (sec)') 

% Stack the figures such that they are not on top of each other
x_wide = .4;
y_high = .4;
x_start = .03;
y_loc = .05;
x_step = .07;

for i = 1:length(fh)
  x_loc = x_start + (i-1) * x_step;
  pos = [x_loc y_loc x_wide y_high];
  set(fh(i),'Units','normalized');
  set(fh(i),'Position',pos);
end % for i = 1:length(fh)

         
