function iono_delay = ionodlay(t_gps,lla,az,el,iono_params)

% iono_delay = ionodlay(t_gps,lla,az,el,iono_params);
%
% Computes the ionospheric group delay for the L1 frequency.  
% This group delay is an apparaent bias in the pseudo-range
% measurements collected by receivers that observe the GPS satellites
% through the ionosphere. This implements the standard GPS ionosphere model
% using the data broadcast in the subframe 4 data (see ICD-GPS 200 for details).
%
% Input:
%   t_gps    - GPS time vector for az/el pairs [GPS_week, GPS_sec] (nx2)
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               valid GPS_sec values are 0-604799
%   lla      - matrix of geodetic position (1x3 or nx3) [lat, lon, height] or
%               (1x2 or nx2) [lat, lon]
%               lat and lon are in rad, height is is not used with 
%               the Klobuchar (GPS-ICD-200) ionosphere model
%               valid latitudes are -pi/2 -> pi/2
%               valid longitudes are -pi -> 2*pi
%   az       - azimuth (-2*pi to 2*pi rad) (nx1)
%   el       - elevation (0 to pi/2 rad) (nx1)
%   iono_params - input data for the ionosphere model (optional)
%                 1x8 matrix in the form 
%                 [alpha_0 alpha_1 alpha_2 alpha_3 beta_0 beta_1 beta_2 beta_3]
%                 where the alpha and beta are transmitted by the GPS satellites.
%                 These parameters are also available from RINEX2 data files
%                 avialable on the World Wide Web at
%                     http://www.ngs.noaa.gov/~don/Data4.html
%                 in the RINEX navigation files (these end with an n).
%                 Default values are 
%                    [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
%                     0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06] 
%                 These default values may or may not be representative of the 
%                 current state of the ionopshere.  Donwload new alpha and beta
%                 values for current times to reflect the current ionosphere.
% Output:
%   iono_dlay   - ionospheric delay (m) (nx1)
%
% See also TROPDLAY 

% Written by: Jimmy LaMance 9/8/97
% Copyright (c) 1998 by Constell, Inc.

% Parkinson et. al, "Global Positioning System: Theory and Applications",
% Vol. 1, pp. 144-149.  There is a typo in the equation for I = T_iono
% in the x^4/4 term.  It should be x^4/24.  See ICD-GPS-200 for correct
% formualtion.

% functions called: ERR_CHK

% WGS-84 constants
RE = 6378137;                  % WGS-84 value in meters
LIGHT_SPEED = 299792458;       % WGS-84 value in m / s

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
iono_delay=[]; 

% Check the number of input arguments and issues a message if invalid
msg = nargchk(4,5,nargin);
if ~isempty(msg)
  fprintf('%s  See help on IONODLAY for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Set the ionosphere model parameters to the default if not provided
if nargin < 5
  iono_params = [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
                 0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06]; 
end % if nargin < 5

estruct.func_name = 'IONODLAY';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 't_gps';
estruct.variable(1).req_dim = [901 2];
estruct.variable(1).var = t_gps;
estruct.variable(1).type = 'GPS_TIME';
  
estruct.variable(2).name = 'lla';
estruct.variable(2).req_dim = [901 3; 901 2; 1 3; 1 2];
estruct.variable(2).var = lla;
estruct.variable(2).type = 'LLA_RAD';
  
estruct.variable(3).name = 'az';
estruct.variable(3).req_dim = [901 1];
estruct.variable(3).var = az;
estruct.variable(3).type = 'ANGLE_RAD';
  
estruct.variable(4).name = 'el';
estruct.variable(4).req_dim = [901 1];
estruct.variable(4).var = el;
estruct.variable(4).type = 'ELEVATION_RAD_0';
  
estruct.variable(5).name = 'iono_params';
estruct.variable(5).req_dim = [1 8];
estruct.variable(5).var = iono_params;
  
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% rename the ionosphere parameters
alpha = iono_params(1:4);
beta = iono_params(5:8);

% allocate the output data matrices and initialize to zero.  This means
% that invalid inputs to the iono model will return no ionosphere delay.                                                                              
iono_delay = zeros(size(az,1),1);

% convert all longitudes to positive (East longitude)
I_az = find(az < 0);
if ~isempty(I_az)
  az(I_az) = az(I_az) + 2 * pi;
end % if ~isempty(I_az)           

% convert the elevation, azimuth, latitude and longitude to from 
% rad to semicircles and cull out the bad values 
semi2rad = pi;
rad2semi = 1 / semi2rad;

el = el * rad2semi;        % elevation angle (semi-circles)
az = az * rad2semi;        % azimuth angle (semi-circles)

lla_temp(:,1:2) = lla(:,1:2) * rad2semi;  % user lat/lon (semi-circles)

t_gps_temp(:,:) = t_gps;

% compute the Earth central angle between the user position and the Earth
% projection of ionosphere intersection point in semicircles (psi)
psi = (0.0137 ./ (el + 0.11)) - 0.022;        % semicircles

% compute the geodetic longitude of the Earth projection of the ionosphere
% intersection point in semicircles
lat_int = lla_temp(:,1) + psi .* cos(az * semi2rad);
I_lat_big = find(lat_int > 0.416);
I_lat_small = find(lat_int < -0.416);

if ~isempty(I_lat_big)
  lat_int(I_lat_big) = ones(length(I_lat_big),1) * 0.416; 
  clear I_lat_big
end % if ~isempty(I_lat_big)

if ~isempty(I_lat_small)
  lat_int(I_lat_small) = -ones(length(I_lat_small),1) * 0.416; 
  clear I_lat_small
end % if ~isempty(I_lat_small)

% compute the geodetic longitude of the Earth projection of the ionosphere
% intersection point in semicircles
lon_int = lla_temp(:,2) + (psi .* sin(az * semi2rad) ./ cos(lat_int * semi2rad));

% compute the geomagnetic latitude of the Earth projection of the ionospheric
% intersection point, mean ionospheric height assumed 350 km, semicircles
phi_m = lat_int + 0.064 * cos((lon_int - 1.617) * semi2rad);

% compute local solar time, t
t = 4.32 * 10^4 * lon_int + t_gps_temp(:,2);

% verify that t is mod 86400 seconds (1 day)
t = rem(t,86400);

% make sure that the t value is positive
I_tneg = find(t < 0);

if any(I_tneg)
  t(I_tneg) = t(I_tneg) + 86400;
end 

if any(find(t < 0))
  fprintf('Found local times less than zero in ionodlay.\n');
  keyboard
end 

% compute the obliquity factor, F
F = 1.0 + 16.0 * (0.53 - el).^3; 

% compute the period, seconds
PER = beta(1) + beta(2) * phi_m + beta(3) * phi_m.^2 + beta(4) * phi_m.^3;

% verify that all the periods (PER) are greater than 72,000 seconds
I_small_period = find(PER < 72000);

% set all period less than 72,000 seconds to 72,000
if ~isempty(I_small_period)
  PER(I_small_period) = 72000 * ones(length(I_small_period),1);
end % if ~isempty(I_small_period)                  

% compute the phase, x (radians)
x = 2 * pi * (t - 50400) ./ PER;

% compute the amplitude
AMP = alpha(1) + alpha(2) * phi_m + alpha(3) * phi_m.^2 + alpha(4) * phi_m.^3;

% verify that the amplitude values (AMP) are greater than zero, if not
% set them to zero
I_neg = find(AMP < 0);

if any(I_neg)
  AMP(I_neg) = zeros(length(I_neg),1);
end % if any(I_neg)

% compute the ionospheric delay in seconds using different formualtions for 
% day and night (day and night are obtained from the phase)
I_day = find(abs(x) < 1.57);
I_night = find(abs(x) >= 1.57);

% if there are daytime observations, do the work
if ~isempty(I_day)
  i_delay(I_day) = F(I_day) .* ...
                   (5.0 * 10e-9 + AMP(I_day) .* ...
                   (1 - x(I_day).^2 ./ 2 + x(I_day).^4 ./ 24));
% Note: This equation is incorrect in "GPS: Theory and Applications".
%       The x^4/24 term has as typo as x^4 / 4.

end % if exist(I_day)                

% if there are nighttime observations, do the work
if ~isempty(I_night)
  i_delay(I_night) = F(I_night) * 5.0e-9;
end % if exist(I_night)                

% convert the delay to meters
iono_delay = i_delay' * LIGHT_SPEED;
                                                      
%%%%% END ALGORITHM CODE %%%%%

% end of IONODLAY
  
