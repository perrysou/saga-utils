function [x_ecef] = lla2ecef(lla)

% [x_ecef] = lla2ecef(lla);
%
% Function to compute ECEF coordinates from 
% geodetic latitude, longitude, and height.
%
% Input:
%   lla    - matrix of geodetic (nx3) [lat, lon, height]
%             lat and lon are in rad, 
%             height is meters above the WGS-84 ellipsoid
%             valid latitudes are -pi/2 -> pi/2
%             valid longitudes are -pi -> 2*pi
% Output:
%   x_ecef - ECEF position in meters (nx3) [x_ecef, y_ecef, z_ecef]
%
% See also ECEF2LLA, ECEF2NED, ECEF2ECI, SIDEREAL 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% Reference: 'GPS: Theory and Practice',
% Hoffman-Wellenhoff, pages 33, 255-257.

% functions called: ERR_CHK

% WGS-84 constants
RE =  6378137;        % equatorial radius of Earth, semi-major axis (m)
f = 1 / 298.2572236;  % WGS-84 flattening
b = 6356752.314;      % ellipsiod semi-minor axis

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_ecef=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on LLA2ECEF for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'LLA2ECEF';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'lla';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = lla;
  estruct.variable(1).type = 'LLA_RAD';
  
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

% allocate the x_ecef output matrix to all inf 
% lla's will be filled in the following code
x_ecef = ones(size(lla)) * inf;

% compute the radius of curvature in the prime vertical
n = RE^2 ./ sqrt(RE^2 * cos(lla(:,1)).^2 + b^2 * sin(lla(:,1)).^2);

% x-component
x_ecef(:,1) = (n + lla(:,3)) .* cos(lla(:,1)) .* cos(lla(:,2)); 

% y-component
x_ecef(:,2) = (n + lla(:,3)) .* cos(lla(:,1)) .* sin(lla(:,2));                    

% z-component
x_ecef(:,3) = ((b^2 / RE^2) * n + lla(:,3)) .* sin(lla(:,1));                      

%%%%% END ALGORITHM CODE %%%%%

% end LLA2ECEF

       
