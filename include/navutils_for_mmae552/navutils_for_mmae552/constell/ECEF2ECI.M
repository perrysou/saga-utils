function [x_eci, v_eci] = ecef2eci(GPS_time, x_ecef, v_ecef)

% For position only conversion
% [x_eci] = ecef2eci(GPS_time, x_ecef);
%                     or
% For position and velocity conversion
% [x_eci, v_eci] = ecef2eci(GPS_time, x_ecef, v_ecef);
%
% Function to convert position and velocity from ECEF to ECI coordinates. 
%
% Input:
%   GPS_time - GPS time (nx2) [GPS_week GPS_sec] 
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               valid GPS_sec values are 0-604799
%   x_ecef   - ECEF position in m (nx3) [x_ecef, y_ecef, z_ecef]
%   v_ecef   - ECEF velocity in m/s (nx3) [vx_ecef, vy_ecef, vz_ecef] (optional)
%   Note: If v_ecef is not provided, position only is converted from 
%         ECEF to ECI. v_eci will return filled with inf if no v_ecef is given.
% Output:
%   x_eci    - ECI position in m (nx3) [x_eci, y_eci, z_eci]
%   v_eci    - ECI velocity in m/s (nx3) [vx_eci, vy_eci, vz_eci]
%
% See also ECI2ECEF, SIDEREAL, ECEF2LLA, ECEF2NED 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPS2UTC, SIDEREAL

% WGS-84 constants
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_eci=[]; v_eci=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECEF2ECI for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 3
  v_ecef = zeros(size(x_ecef,1),3);
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ECEF2ECI';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'GPS_time';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = GPS_time;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'x_ecef';
  estruct.variable(2).req_dim = [901 3];
  estruct.variable(2).var = x_ecef;
  
  estruct.variable(3).name = 'v_ecef';
  estruct.variable(3).req_dim = [901 3];
  estruct.variable(3).var = v_ecef;
  
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

% allocate the output x_eci matrix, filled with infinity
% valid inputs will be filled correctly in the following code
x_eci = ones(size(x_ecef)) * inf;
v_eci = ones(size(x_ecef)) * inf;

% compute the current Greenwich sidereal time from the current utc_time

utc_time = gps2utc(GPS_time);    % convert to UTC for sidereal calcs 
omega = sidereal(utc_time);                     % omega is the Sidereal time (rad) 

% now do the coordinate transformation one column at a time (vectorized)
x_eci(:,1) = x_ecef(:,1) .* cos(omega) - x_ecef(:,2) .* sin(omega);
x_eci(:,2) = x_ecef(:,1) .* sin(omega) + x_ecef(:,2) .* cos(omega);
x_eci(:,3) = x_ecef(:,3);

if nargin == 3   % do velocity rotation also

  % first convert ECEF velocity to ECI frame
  % these are ECEF velocities written in the ECI frame
  v_eci(:,1) = v_ecef(:,1) .* cos(omega) - v_ecef(:,2) .* sin(omega);
  v_eci(:,2) = v_ecef(:,1) .* sin(omega) + v_ecef(:,2) .* cos(omega);
  v_eci(:,3) = v_ecef(:,3);

  % compute omega cross r term to convert from ECEF to inertial velocities
  e_rate_eci = [0 0 EARTH_RATE]' * ones(1,size(x_eci(:,:),1));
  o_cross_r = cross(e_rate_eci,x_eci(:,:)')';

  % convert from ECEF to ECI velocity by adding the omega cross r term
  % (these calculations are done in the eci frame)
  v_eci(:,:) = v_eci(:,:) + o_cross_r;

end % if nargin == 3        

%%%%% END ALGORITHM CODE %%%%%

% end ECEF2ECI

