function [x_ecef, v_ecef] = eci2ecef(GPS_time, x_eci, v_eci)

% Function to convert position and velocity from ECI to ECEF coordinates. 
%
% For position only conversion
% [x_ecef] = eci2ecef(GPS_time, x_eci);
%                     or
% For position and velocity conversion
% [x_ecef, v_ecef] = eci2ecef(GPS_time, x_eci, v_eci);
%
% Input:
%   GPS_time - GPS time (nx2) [GPS_week GPS_sec]
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               valid GPS_sec values are 0-604799
%   x_eci    - ECI position in meters (nx3) [xi, yi, zi]
%   v_eci    - ECI velocity in m/s (nx3) [vxi, vyi, vzi] (optional)
%   Note: If v_eci is not provided, position only is converted from 
%         ECI to ECEF. v_ecef will return filled with inf if no v_eci is given.
% Output:
%   x_ecef   - ECEF position in meters (nx3) [xe, ye, ze]
%   v_ecef   - ECEF velocity in m/s (nx3) [vxe, vye, vze] 
%
% See also ECEF2ECI, SIDEREAL, ECEF2LLA, ECEF2NED 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPS2UTC, SIDEREAL

% WGS-84 constants
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_ecef=[]; v_ecef=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECI2ECEF for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 3
  v_eci = zeros(size(x_eci,1),3);
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ECI2ECEF';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'GPS_time';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = GPS_time;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'x_eci';
  estruct.variable(2).req_dim = [901 3];
  estruct.variable(2).var = x_eci;
  
  estruct.variable(3).name = 'v_eci';
  estruct.variable(3).req_dim = [901 3];
  estruct.variable(3).var = v_eci;
  
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

% allocate the output x_ecef matrix, filled with infintiy
% valid inputs will be filled correctly in the following code
x_ecef = ones(size(x_eci)) * inf;
v_ecef = ones(size(x_eci)) * inf;

% compute the current Greenwich sidereal time from the current utc_time

utc_time = gps2utc(GPS_time(:,:));    % convert to UTC for sidereal calcs
omega = sidereal(utc_time);                     % omega is the Sidereal time (rad) 

% convert from ECI to ECEF
x_ecef(:,1) = x_eci(:,1) .* cos(omega) + x_eci(:,2) .* sin(omega);
x_ecef(:,2) = -x_eci(:,1) .* sin(omega) + x_eci(:,2) .* cos(omega);
x_ecef(:,3) = x_eci(:,3);

if nargin == 3   % do velocity rotation also

  % First, subtract Earth rate term in ECI frame to get ECEF velocity
  % in ECI frame.
  % Second, convert ECEF velocity from ECI to ECEF frame

  e_rate_eci = [0 0 EARTH_RATE]' * ones(1,size(x_eci(:,:),1));
  o_cross_r = cross(e_rate_eci,x_eci(:,:)')';

  % convert to ECEF velocity by subtracting the omega cross r term
  % (do a software no-no and use the same variable name for two different
  %  physical quantities to save storage space)

  v_eci(:,:) = v_eci(:,:) - o_cross_r;  % v_eci is now v_ecef in i

  % now convert ECEF velocity from ECI to ECEF frame
  v_ecef(:,1) = v_eci(:,1) .* cos(omega) + v_eci(:,2) .* sin(omega);
  v_ecef(:,2) = -v_eci(:,1) .* sin(omega) + v_eci(:,2) .* cos(omega);
  v_ecef(:,3) = v_eci(:,3);

end % if nargin == 3        

%%%%% END ALGORITHM CODE %%%%%

% end of ECI2ECEF
