function [t_dop,prn,doppler,dop_orb,dop_err]=...
       doppler(t_user,x_user,v_user,t_gps,x_gps,v_gps,model,seed,dop_noise);

% [t_dop, prn, doppler] = doppler(t_user, x_user, v_user, t_gps, x_gps, v_gps);
%                          or
% [t_dop,prn,doppler,dop_orb,dop_err]=...
%      doppler(t_user,x_user,v_user,t_gps,x_gps,v_gps,model,seed,dop_noise);
%
% Computes Doppler measurements given a user trajectory and the GPS/GLONASS 
% satellite positions and velocities. Errors that can be included in 
% the Doppler measurement include S/A (dither), receiver clock bias and drift, 
% and receiver Doppler measurement noise.
%
% Input:
%   t_user - GPS time vector for user trajectory [GPS_week, GPS_sec] (nx2)
%             valid GPS_week values are 1-3640 (years 1980-2050)
%             valid GPS_sec values are 0-604799
%   x_user - ECEF/ECI position vectors for the user vehicle [x,y,z] (nx3) (m)
%   v_user - ECEF/ECI vel. vectors for the user vehicle [xv,yv,zv] (nx3) (m/s)
%   t_gps  - GPS time vector for GPS positions [GPS_week, GPS_sec] (mx2)
%             valid GPS_week values are 1-3640 (years 1980-2050)
%             valid GPS_sec values are 0-604799
%   x_gps  - ECEF/ECI position vectors for GPS satellites [prn,x,y,z] (mx4) (m)
%   v_gps  - ECEF/ECI velocity vectors for GPS satellites [xv,yv,zv] (mx3) (m/s)
%   model  - flags controlling which contributions to the Doppler errors are 
%             modeled (optional).  (1x3) [sa_dither user_clock receiver_noise]
%             a value of 1 indicates useage of the model and a value of zero 
%             indicates no use of that model. Default = [1 1 1].
%   seed   - seed value for random number generator (optional). Default = 0. 
%   dop_noise - 1 sigma estimate of the receiver Doppler noise (optional) (1x1)
%                (m/s). Default = .3 m/s.
% Output:
%   t_dop   - GPS time associated with the Doppler measurement, 
%              [GPS_week GPS_sec] (kx2), 
%              k = num_time_steps x number of visible satellites
%   prn     - satellite number for this Doppler measurement (kx1)
%   doppler - Doppler measurements associated with the corresponding 
%              t_dop and prn (kx1) (m/s)
%   dop_orb - GPS/GLONASS satellite orbits associated with this 
%              measurement (kx6) [orb_x orb_y orb_z orb_vx orb_vy orb_vz]
%   dop_err - modeled errors added to the Doppler measurement (m/s) (kx3)
%                [sa_dither user_clock receiver_noise]
%
%   Note: Doppler measurements will be generated at valid t_gps times.  The user
%         trajectory will be interpolated to the t_gps times if they are not
%         coincident with the t_gps times on input.
%
% See also PSEUDO_R, SA_EPS, SA_CLOCK, CLOCKERR, TROPDLAY, IONODLAY

% Written by: Jimmy LaMance 9/2/97
% Copyright (c) 1998 by Constell, Inc.

% Reference: 'GPS: Theory and Practice',
%             Hoffman-Wellenhoff, pages 92-93, 182. 
%
%             "Global Positioning System: Theory and Applications", 
%             Volume 1, Parkinson and Spilker, pages 411-412.

% functions called: SA_CLOCK, CLOCKERR, NORMVECT, ERR_CHK, LOS

% WGS-84 constants
LIGHT_SPEED = 299792458;       % WGS-84 value in m / s
EARTH_RATE = 7.2921151467e-5;  % WGS-84 value in rad / s 

% GPS constants
L1_FREQ = 1575.42e6;                      % Hz (1575.42 MHz)
L1_WAVELENGTH = LIGHT_SPEED / L1_FREQ;

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_dop=[]; prn=[]; doppler=[]; dop_orb=[]; dop_err=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(6,9,nargin);
if ~isempty(msg)
  fprintf('%s  See help on DOPPLER for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 7
  model = [1 1 1];     % set to the default value
end % if nargin < 7 

if nargin < 8
  seed = 0;     % set to the default value
end % if nargin < 8

% verify that dop_noise is a 1x1, if provided, if not set to the default
if nargin < 9 
  dop_noise = .3;     % set to the default value 0.3 m/s
end % if nargin < 9 

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'DOPPLER';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_user';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_user;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'x_user';
  estruct.variable(2).req_dim = [901 3];
  estruct.variable(2).var = x_user;
  
  estruct.variable(3).name = 'v_user';
  estruct.variable(3).req_dim = [901 3];
  estruct.variable(3).var = v_user;
  
  estruct.variable(4).name = 't_gps';
  estruct.variable(4).req_dim = [902 2];
  estruct.variable(4).var = t_gps;
  estruct.variable(4).type = 'GPS_TIME';
  
  estruct.variable(5).name = 'x_gps';
  estruct.variable(5).req_dim = [902 4];
  estruct.variable(5).var = x_gps;

  estruct.variable(6).name = 'v_gps';
  estruct.variable(6).req_dim = [902 3];
  estruct.variable(6).var = v_gps;

  estruct.variable(7).name = 'model';
  estruct.variable(7).req_dim = [1 3];
  estruct.variable(7).var = model;
  
  estruct.variable(8).name = 'seed';
  estruct.variable(8).req_dim = [1 1];
  estruct.variable(8).var = seed;
  
  estruct.variable(9).name = 'dop_noise';
  estruct.variable(9).req_dim = [1 1];
  estruct.variable(9).var = dop_noise;
  
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
  
% Compute geometric los vectors (range vector)
% Ths indices will be used to sort the remainder of the data
[t_dop,los,los_ndex] = los(t_user,x_user,t_gps,x_gps);

[los_norm, pr] = normvect(los); 

I_gps = los_ndex(:,2);
I_user = los_ndex(:,1);
prn = x_gps(I_gps,1); 
dop_orb = [x_gps(I_gps,2:4) v_gps(I_gps,:)];

% compute the relative velocity of the user trajectory and the GPS satellite
vel = v_gps(I_gps,:) - v_user(I_user,:);

% compute the Doppler shift as the dot product of the LOS with the velocity
doppler = dot(los_norm',vel')';

% add dither effect
if model(1) == 1 
  % set the sigma value for the SA clock model to the RTCM default values
  sa_model_data = [23 .28 118];    % sigma PR, sigma PR-rate, tau (decorrelation
                                   % time) (m, m/s, dimensionless)                                          

  % compute the dither effects using the default model (RTCA parameters)
  [sa_clk_err, sa_clk_drift] = sa_clock(t_dop, prn, sa_model_data, seed);
  doppler = doppler + sa_clk_drift;                 
  
else
  sa_clk_drift = zeros(size(doppler,1),1);  
end % model(1) == 1  

% compute user clock bias and drift
if model(2) == 1
  % set the value for the clock model to the default
  clk_model = [4e-19 1.58e-18];     % bias and frequency noise values
  
  % compute the clock bias using the default model 
  % roughly corresponding to a crystal oscillator
  [clk_bias clk_drift] = clockerr(t_dop, clk_model, seed);

  % add the clock bias to the pseudo-range and accumulated phase measurements
  doppler = doppler + clk_drift;                 

else
  clk_drift = zeros(size(pr,1),1);  
end % if model(2) == 1

% add receiver white noise 
if model(3) == 1
  dop_white = randn(size(doppler,1),1) * dop_noise;

  doppler = doppler + dop_white; 
   
else
  dop_white = zeros(size(doppler,1),1);
  
end % if model(6) == 1

% build up the PR error matrix for output
dop_err = [sa_clk_drift clk_drift dop_white];

%%%%% END ALGORITHM CODE %%%%%

% end of DOPPLER
