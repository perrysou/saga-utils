% ex_dgps.m
%
% Example script for generating GPS measurements and computing position
% and velocity solutions using the raw data. This example is for a static 
% receiver fixed to the surface of the Earth with a base station located
% approximately 250 miles away.

% Written by: Jimmy LaMance 
% Copyright (c) 1998 by Constell, Inc.

% functions called: UTC2GPS, FIND_ALM, READYUMA, ALM2GEPH,
%                   PROPGEPH, ECEF2NED, PSEUDO_R, ECEF2LLA,
%                   NED2AZEL, VIS_DATA, DOPPLER, LSNAV, DIFFCORR,
%                   ADD_DPR, LOS, NED2DOPS, DOPS2ERR, PLOTPASS

% Variables to be used in the examples
start_time = [1998 7 4 0 0 0];        % UTC start time
duration = 1*3600;                     % 1 hour duration
base_station = [40*pi/180 255*pi/180 2000];       % base station coordinates
remote_station = [42*pi/180 254*pi/180 2000];     % remote station coordinates

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
t_remote = t_base;

% Convert the base and remote station coordinates to ECEF
x_base = ones(size(t_base,1),1) * lla2ecef(base_station);  
v_base = zeros(size(x_base));
x_remote = ones(size(t_remote,1),1) * lla2ecef(remote_station);
v_remote = zeros(size(x_remote));

% Compute the base line length
base_line_ecef = x_base(1,:) - x_remote(1,:);
base_line_ned = ecef2ned(base_line_ecef,ecef2lla(x_base(1,:)));
base_line_km = norm(base_line_ned) / 1000;

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
pr_errors_base = pr_errors_base(I_pass,:);
base_ndex = base_ndex(I_pass,:);

% Rename the pr_errors for easier use later
pr_sa_errb = pr_errors_base(:,1) + pr_errors_base(:,2);
trop_base = pr_errors_base(:,3) + pr_errors_base(:,4);
iono_base = pr_errors_base(:,5);
clk_biasb = pr_errors_base(:,6);
clk_driftb = pr_errors_base(:,7);

remote_mask = 5*d2r;           
% Model SA, tropo, iono, and receiver clock.  Do not model satellite motion,
% Earth rotation, satellite clocks, line baises, or relativity.  The LSNAV
% function used to compute navigation solutions does not currently support
% this PR modeling.  However, these effects could be modeled and then absorbed
% into the differential correction.
remote_model = [1 1 1 1 1 1 0 0 0 0 0];
remote_code_noise = 1;
remote_carrier_noise = .01;
remote_seed = 0;
[t_pr_remote,prn_remote,pr_remote,pr_errors_remote,remote_ndex] = ...
               pseudo_r(t_remote,x_remote,v_remote,t_gps,[prn_gps x_gps],v_gps,...
                        remote_model,remote_seed,...
                        remote_code_noise, remote_carrier_noise,gps_ephem);

% Save the remote orbits associate with this PR data into it's own
% matrix for easier booking keeping later.
pr_orb_remote = x_gps(remote_ndex(:,2),:);

% Compute Doppler measurements at the remote site to be used 
% in computing velocity solutions
remote_doppler_model = [1 1 1]; 
remote_dop_noise = .3;
[t_dop,prn_remote_dop,dopp,dop_orb,dop_err] = ...
     doppler(t_remote,x_remote,v_remote,t_gps,[prn_gps x_gps],v_gps,...
             remote_doppler_model,remote_seed,remote_dop_noise);

% Compute masking for the remote receiver.  Here again, the NED system is used.
los_remote_ecef = x_gps(remote_ndex(:,2),:) - x_remote(remote_ndex(:,1),:);

% Rotate the LOS to NED, using the base station as the reference for the NED
% coordinate system and rotation.
ref_lla = ecef2lla(x_remote(remote_ndex(:,1),:));

los_base_ned = ecef2ned(los_base_ecef,ref_lla);

% Compute the azimuth/elevation of the NED vectors
[az, el] = ned2azel(los_base_ned);

% Apply the masking in the NED coordinate system
[visible_data, I_pass] = vis_data(remote_mask,[az el]);

% Remove all of the base station data that did not pass the masking test
t_pr_remote = t_pr_remote(I_pass,:);
prn_remote = prn_remote(I_pass,:);
pr_remote = pr_remote(I_pass,:);
pr_orb_remote = pr_orb_remote(I_pass,:);
pr_errors_remote = pr_errors_remote(I_pass,:);
remote_ndex = remote_ndex(I_pass,:);
dopp = dopp(I_pass,:);
dop_orb = dop_orb(I_pass,:);

% Rename the pr_errors for easier use later
pr_sa_errr = pr_errors_remote(:,1) + pr_errors_remote(:,2);
trop_remote = pr_errors_remote(:,3) + pr_errors_remote(:,4);
iono_remote = pr_errors_remote(:,5);
clk_biasr = pr_errors_remote(:,6);
clk_driftr = pr_errors_remote(:,7);

% Compute a position and velocity solution at the remote site without DGPS
[t_nav,x_nav,num_remote,nav_index,v_nav] = ...
     lsnav(t_pr_remote,pr_remote(:,1),...
           [prn_remote x_gps(remote_ndex(:,2),:)],[x_remote(1,:) 0],dopp, ...
            dop_orb(:,4:6),[0 0 0 0]); 

% Compute differential corrections given the base station locations.  Model
% the GPS satellite clocks, relativistic effects, satellite motion, and
% Earth rotation.  SA, ionosphere, troposphere, and receiver clock are not
% modeled in the PR for computing differential corrections.  If a term is 
% not modeled in the correction, but was included in the PR measurement, it
% will effectively be absorbed into the differential correction. 
diff_model = [0 0 0 0 0 0 0 0 0 0 0];   
[t_dpr, dpr] = ...
       diffcorr(t_pr_base, [prn_base pr_base], x_gps(base_ndex(:,2),:), ...
                v_gps(base_ndex(:,2),:), x_base(1,:),diff_model,gps_ephem);
       
% Set the differential correction rate of change to zero
dprr = zeros(size(dpr,1),1);

% Add the corrections to the remote PR measurements
% the return variable cpr_index is an index into the input remote
% site pseudoranges that have differential corrections applied to them  
[t_cpr, prn_cpr, cpr, cpr_index] = ...
      add_dpr(t_pr_remote,[prn_remote pr_remote],t_dpr,dpr,dprr,max_latency);

% Compute a LS nav solotuion at the remote site with DGPS
[t_nav_dgps,x_nav_dgps,num_dgps] = ...
     lsnav(t_cpr,cpr(:,1),[prn_remote(cpr_index) pr_orb_remote(cpr_index,:)],...
           [x_remote(1,:)+50 0]); 

% Compute residual troposphere errors using the ADD_DPR function.  this function
% handles all of the common visibility analysis.  Instead of pseudoranges and
% pseudorange corrections, we will send it troposphere errors.  To obtain the
% troposphere differences, we will send in the negative of one set of tropo
% errors.
[t_dtrop, prn_dtrop, dtrop, dtrop_index] = ...
      add_dpr(t_pr_remote,[prn_remote trop_remote],t_pr_base,...
              [prn_base -trop_base],dprr,0);

% Compute the residual ionospheric errors the same way
[t_diono, prn_diono, diono, diono_index] = ...
      add_dpr(t_pr_remote,[prn_remote iono_remote],t_pr_base,...
              [prn_base -iono_base],dprr,0);

% Compute LOS at remote station in ECEF frame
[t_los, los, los_ind] = los(t_gps(1,:), x_remote(1,:), t_gps, [prn_gps x_gps]);

veh_num = prn_gps(los_ind(:,2));

% Rotate ECEF los to NED
[x_ned] = ecef2ned(los, ecef2lla(x_remote(1,:)));

% Compute masking
[az, el] = ned2azel(x_ned);
[az_el, I_pass] = vis_data(remote_mask, [az el]);

% Compute DOPs at the remote station (with masking for the remote station)
[remote_dops, t_dops, num_base] = ned2dops(x_ned(I_pass,:),t_los(I_pass,:));

% Estimate position errors from DOPs
sigma_pr = 40;
[pos_err_dops] = dops2err(remote_dops,sigma_pr);
     
% Compute the position errors with and without DGPS
% these vectors are in ECEF
pos_err = x_nav(:,1:3) - ones(size(x_nav,1),1) * x_remote;
pos_err_dgps = x_nav_dgps(:,1:3) - ones(size(x_nav,1),1) * x_remote(:,1:3);

% Rotate position difference to NED relative to the truth location
% for easier interpretation                                       
[pos_err_ned] = ecef2ned(pos_err, ecef2lla(x_remote));
[pos_err_ned_dgps] = ecef2ned(pos_err_dgps, ecef2lla(x_remote));

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
  sprintf('Position Errors without DGPS with a %d km Baseline',...
           round(base_line_km));
title(title_string) 

subplot(2,1,2)
plot(pos_err_ned_dgps)
legend('North','East','Down',-1);
ylabel('Position Err (m)')
title_string = ...
  sprintf('Position Errors with DGPS with a %d km Baseline',...
           round(base_line_km));
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

% plot differential troposphere contribution to the differential PR error
fh(4) = plotpass(t_dtrop,dtrop,prn_dtrop,...
         'Example of Differential Troposphere Effects on Pseudorange Measurements',...
         'Delta Tropo (m)');

% plot ionosphere contribution to the PR error
fh(5) = plotpass(t_pr_base,iono_base,prn_base,...
         'Example of Ionosphere Effects on Pseudorange Measurements',...
         'Iono (m)');

% plot differential ionosphere contribution to the differential PR error
fh(6) = plotpass(t_diono,diono,prn_diono,...
         'Example of Differential Ionosphere Effects on Pseudorange Measurements',...
         'Delta Iono (m)');

% plot the receiver clock bias for a single satellite 
% the receiver clock bias is common to all satellites
I = find(prn_base == prn_base(2));
fh(7) = plotpass(t_pr_base(I,:),clk_biasb(I),prn_base(I),...
         'Example of Receiver Clock Bias on Pseudorange Measurements',...
         'Clock Bias (m)');

fh(8) = plotpass(t_pr_base(I,:),clk_driftb(I),prn_base(I),...
         'Example of Receiver Clock Drift on Pseudorange Measurements',...
         'Clock Drift (m/2)');

% generate a plot of the resulting position error from the DOPs
t_dop_plot = (t_dops(:,1) * 604800 + t_dops(:,2)) - ...
             (t_dops(1,1) * 604800 + t_dops(1,2));

fh(9) = figure;
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
end % for i = 1:legnth(fh)

if exist('save_plot_data')
  % save the data to a mat file for regeneration later
  save demoplt1 plot_time pos_err_ned pos_err_ned_dgps pr_base ...
       t_pr_base pr_sa_errb clk_biasb clk_driftb trop_base trop_remote ...
       iono_base iono_remote ...
       prn_base base_line_km remote_dops t_dop_plot pos_err_dops ...
       t_nav num_remote num_base num_dgps ...
       t_diono prn_diono diono ...
       t_dtrop prn_dtrop dtrop
end % if save_plot_data == 1

         