function [t_out, prn, x, v] = propgeph(gps_ephem, t_start, t_stop, dt)

% [t_out, prn, x, v] = propgeph(gps_ephem, t_start, t_stop, dt);
%
% Computes satellite positions in ECEF coordinates from a GPS ephemeris.  
%
% Input:
%   gps_ephem - ephemeris matirx for all satellites (nx24), with columns of
%                [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                Ephemeris parameters are from ICD-GPS-200 with the 
%                exception of perigee_rate.
%   t_start   - start time in GPS format (1x2) or (nx2) [GPS_week GPS_sec] 
%                If t_stop is not provided, orbits will be computed at each
%                of the t_start times. Must be 1x2 if t_stop is given.
%                valid GPS_week values are 1-3640 (years 1980-2050)
%                valid GPS_sec values are 0-604799
%   t_stop    - stop time in GPS time format (1x2) [GPS_week GPS_sec] (optional),
%                default = 600 seconds past start
%   dt        - output interval, seconds (optional), default = 2 sec
%
% Output:
%   t_out     - GPS time vector output [GPS_week, GPS_sec] (nx2)
%   prn       - S/C (PRN) number (nx1)
%   x         - S/C position in ECEF meters (nx3)
%   v         - S/C velocity in ECEF meters/second (nx3)
%
% Note: Output are in time order with all of the satellites at time #1
%       first, followed by the satellites at time #2, etc.  The position and
%       velocity correspond to the times and satellite numbers in the 
%       corresponding columns. The t_out and prn matrices have the following form
%                           wk  s  prn   X    Y    Z    Vx    Vy    Vz
%        [t_out prn x v] = [933 1   1   x11  y11  z11  vx11  vy11  vz11
%                           933 1   2   x12  y12  z12  vx12  vy12  vz12
%                            .      .    .    .    .     .     .     .
%                            .      .    .    .    .     .     .     .
%                           933 1   24  x124 y124 z124 vx124 vy124 vz124
%                           933 2   1   x21  y21  z21  vx21  vy21  vz21
%                           933 2   2   x22  y22  z22  vx22  vy22  vz22
%                            .      .    .    .    .     .     .     .
%                           933 600 24   .    .    .     .     .     .  ]
%
% See also ALM2GEPH, KEP2GEPH

% Written by: Jimmy LaMance 10/18/96
% Copyright (c) 1998 by Constell, Inc.

% This is taken directly from
% "Global Positioning System: Theory and Applications", Parkinson and Spilker,
% Vol. 1, pages 132-138 with the equation for x(:,2) (yk in Table 8) corrected.  

% functions called: ERR_CHK, KEPLR_EQ

% WGS-84 constants
MU_EARTH = 3.986005e14;       % WGS-84 value in m^3/s^2
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 
PI = 3.1415926535898;         % accepted GPS value for pi

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_out=[]; prn=[]; x=[]; v=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,4,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PROPGEPH for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in optional variable
if nargin < 3 & size(t_start,2) == 2
  duration = 600;
  start_sec = gpst2sec(t_start);
  stop_sec = start_sec+duration;
  t_stop = sec2gpst(stop_sec);             
end

if nargin < 4
  dt = 2;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'PROPGEPH';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'gps_ephem';
  estruct.variable(1).req_dim = [901 24];
  estruct.variable(1).var = gps_ephem;
  
  estruct.variable(2).name = 't_start';
  estruct.variable(2).req_dim = [1 2; 902 2];
  estruct.variable(2).var = t_start;
  estruct.variable(2).type = 'GPS_TIME';
  
  estruct.variable(3).name = 't_stop';
  estruct.variable(3).req_dim = [1 2; 902 2];
  estruct.variable(3).var = t_stop;
  
  estruct.variable(4).name = 'dt';
  estruct.variable(4).req_dim = [1 1];
  estruct.variable(4).var = dt;
  
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
% find the earliest starting week in the ephemeris data and start data
week_min = min([min(gps_ephem(:,19)) t_start(1)]);

% compute the start and stop times
stop_time_seconds = (t_stop(1) - week_min) * 86400 * 7 + t_stop(2); 
start_time_seconds = (t_start(1) - week_min) * 86400 * 7 + t_start(2); 

% build the linear time series
if nargin < 3
  t_input = (t_start(:,1)  - week_min) * 86400 * 7 + t_start(:,2);
  t_input = t_input';       % get the t_input matrix in the right shape
else
  t_input = [start_time_seconds:dt:stop_time_seconds];
end

% Search the prn column and see how many satellites are in the ephemeris.
for i = 1:size(gps_ephem,1)
  prn(i) = gps_ephem(i,1);
end % for i = 1:length(gps_ephem,1)

num_active_satellites = length(prn);     % we know we have at most this many
% search and sort for mutiple satellites with different ephemeris and 
% subtract the multiples from the total number
for i = 1:length(prn)
  I_temp = find(gps_ephem(:,1) == prn(i));
  num_active_satellite = num_active_satellites - (length(I_temp) - 1);
  clear I_temp
end % for

% loop over the number of active satellites just found
% find how many ephemerides are given for each satellite (possible to have 
% multiple ephemeris for 1 satellite)
for i = 1:num_active_satellite
  num_prn(i) = length(find(gps_ephem(:,1) == prn(i)));                                                                                                                 % different ephemeris for different times 
end % for i = 1:length(prn)                                                  

% find the number of time steps to be evaluated
num_time_steps = size(t_input,2);  

% get total size for memory allocation
total_size = num_time_steps * num_active_satellites; 

% allocate output matrix space
x = ones(total_size,3) * NaN;
v = ones(total_size,3) * NaN;
o_cross_r = ones(total_size,3) * NaN;
t_out = ones(total_size,2) * NaN;
prn = ones(total_size,1) * NaN; 

for jjj = 1:num_active_satellites
  output_index = [jjj:num_active_satellites:total_size]; % indices for output
  
  % convert the ephemeris ephoch into 
  % the same time scale used elsewhere (sec past min. week)
  ephemeris_epoch = (gps_ephem(jjj,19) - week_min) * 86400 * 7 + ...
                     gps_ephem(jjj,17);
  t = t_input - ephemeris_epoch; % time past ephemeris epoch
                                 % Toe = gps_ephem(:,17)
  a = gps_ephem(jjj,5).^2;       % compute a, sqrt(a) = gps_ephem(:,5) (m)
  e = gps_ephem(jjj,4);          % load eccentricity to a local variable
  n0 = sqrt(MU_EARTH / a.^3);    % compute the nominal mean motion (rad/s)
  n = n0 + gps_ephem(jjj,3);     % corrected mean montion
                                 % delta_n = gps_ephem(:,3) 
  M = gps_ephem(jjj,2) + n .* t; % compute mean anomaly, M0 = gps_ephem(:,2)   

  % load perigee to a local variable and add perigee rate
  perigee  = gps_ephem(jjj,8) + gps_ephem(jjj,24) .* t;      

  % check input values of e to verify it is an eccentric orbit
  clear I
  I = find(e > .6);
  if any(I)
    fprintf('Warning message from propgeph.m ...\n')
    fprintf('Results may be in error for eccentricities greater than 0.6.\n')
    fprintf('This propagator is not designed for parabolic or hyperbolic orbits.\n')
    fprintf('Results may be invalid.  e = %f.\n \n',e(I) )
  end % if
  
  % solve Keplers eq for eccentric anomaly (E) from M = E - e sin(E)
  [E] = keplr_eq(M,e);
  
  % compute true anomaly
  nu = atan2( sqrt(1 - e.^2) .* sin(E) ./ (1 - e .* cos(E)), ...
              (cos(E) - e) ./ (1 - e .* cos(E)) ); 

  % compute the argument of latitude
  phi = nu + perigee;

  % compute corrections to argument of latitude

  % compute sinusoidal terms from the ephemeris data
  del_u = gps_ephem(jjj,12) .* sin(2 .* phi) + ... % Cus = gps_ephem(jjj,12)
          gps_ephem(jjj,11) .* cos(2 .* phi);      % Cuc = gps_ephem(jjj,11)

              
  del_r = gps_ephem(jjj,14) .* sin(2 .* phi) + ... % Crs = gps_ephem(jjj,14)
          gps_ephem(jjj,13) .* cos(2 .* phi);      % Crc = gps_ephem(jjj,13)
           
  del_i = gps_ephem(jjj,16) .* sin(2 .* phi) + ... % Cis = gps_ephem(jjj,16)
          gps_ephem(jjj,15) .* cos(2 .* phi);      % Cic = gps_ephem(jjj,15)

  u = phi + del_u;                      % argument of latitude correction
  r = a .* (1 - e .* cos(E)) + del_r;   % radius correction 

  % inclination correction
  %i = gps_ephem(jjj,7), i_rate = gps_ephem(jjj,10)
  i = gps_ephem(jjj,7) + del_i + gps_ephem(jjj,10) .* t; 

  xo = (r .* cos(u))';    % satellite x-position in orbital plane
  yo = (r .* sin(u))';    % satellite y-position in orbital plane

  % corrected longitude of ascending node for node rate and Earth rotation
  % Asc node = gps_ephem(jjj,6), node rate = gps_ephem(jjj,9)
  node = gps_ephem(jjj,6) + (gps_ephem(jjj,9) - EARTH_RATE) .* t -  ...
         EARTH_RATE * gps_ephem(jjj,17);   %   Toe = gps_ephem(jjj,17)

  % satellite x-position in ECEF (m)
  x(output_index,1) = xo .* cos(node') - yo .* cos(i') .* sin(node');    
  % satellite y-position in ECEF (m)
  x(output_index,2) = xo .* sin(node') + yo .* cos(i') .* cos(node');    
  % satellite z-position in ECEF (m)
  x(output_index,3) = yo .* sin(i');                                     

  % compute satellite velocity in the perifocal coordinate system, 
  % from Bate, Mueller, and White pages 73 & 82
  r = sqrt(x(output_index,1).^2 + x(output_index,2).^2 + x(output_index,3).^2);
  
  p = r .* (1 + e .* cos(nu'));

  vp = -sqrt(MU_EARTH ./ p) .* sin(nu');
  vq = sqrt(MU_EARTH ./ p) .* (e + cos(nu'));

  % inertial velocity vector in ECEF coordinates
  v(output_index,1) = vp .* (cos(node') .* cos(perigee') - ...
                      sin(node') .* sin(perigee') .* cos(i')) + ...
                      vq .* (-cos(node') .* sin(perigee') - ...
                      sin(node') .* cos(perigee') .* cos(i')); 
  v(output_index,2) = vp .* (sin(node') .* cos(perigee') + ...
                      cos(node') .* sin(perigee') .* cos(i')) + ...
                      vq .* (-sin(node') .* sin(perigee') + ...
                      cos(node') .* cos(perigee') .* cos(i'));
  v(output_index,3) = vp .* sin(perigee') .* sin(i') + ...
                      vq .* cos(perigee') .* sin(i');          
  
  prn(output_index) = ones(size(output_index)) * gps_ephem(jjj,1);
  t_out(output_index,2) = t_input';
                                                                                       
end % for

if ~any(t_out)
  fprintf('Error message from propgeph.m\n')
  fprintf('No active satellite found.\n')
  fprintf('Possible invalid return arguments. Check inputs.\n\n')
  return
end % if

% compute omega cross r term for use in converting from ECEF to inertial
% velocities
e_rate_ecef = [0 0 EARTH_RATE]' * ones(1,size(x,1));
o_cross_r = cross(e_rate_ecef,x')';

% convert to ECEF velocity by removing the omega cross r term
v = v - o_cross_r;        

% reformat the time matrix (get it back into GPS time)

% fill in the time matrix with this week value
t_out(:,1) = ones(size(t_out,1),1) * week_min; 

I_roll = find(t_out(:,2) >= 86400 * 7);  % check for and find week roll-overs
                                                   
if any(I_roll)   % if there is(are) a roll-over(s), convert it to weeks
  % fix the week
  t_out(I_roll,1) = t_out(I_roll,1) + fix(t_out(I_roll,2) ./ (86400 * 7));    
  % fix the seconds
  t_out(I_roll,2) = rem(t_out(I_roll,2), (86400 * 7));                        
end % if any(I_roll)
  
%%%%% END ALGORITHM CODE %%%%%

% end PROPGEPH      
