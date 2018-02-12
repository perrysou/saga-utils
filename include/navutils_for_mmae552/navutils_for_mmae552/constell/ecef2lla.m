function [lla] = ecef2lla(x_ecef, method)

% [lla] = ecef2lla(x_ecef, method);
%
% Function to compute geodetic latitude, longitude, and 
% height from ECEF coordinates.
%
% Input:
%   x_ecef - ECEF position in meters (nx3) [x_ecef, y_ecef, z_ecef]
%   method - flag indicating method to be used, optional
%             0 - fast but approx. (lat, long errors < 1e-4 rad, alt_err < 100m)
%                       (default) 
%             1 - slower (iterative method, more accurate)
%
% Output:
%   lla    - geodetic lat, lon, and height (nx3) [lat, lon, height]
%             lat and lon are in rad (mod 2*pi), 
%             height is meters above the WGS-84 ellipsoid
%
% See also LLA2ECEF, ECEF2NED, ECEF2ECI, SIDEREAL 

% Written by: Jimmy LaMance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% Reference: 'GPS: Theory and Practice',
% Hoffman-Wellenhoff, pages 33, 255-257.

% functions called: ERR_CHK

% WGS-84 constants
RE =  6378137;          % equatorial radius of Earth, semi-major axis (m)
f = 1 / 298.257223563;  % WGS-84 flattening
b = 6356752.314;        % ellipsiod semi-minor axis

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
lla=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECEF2LLA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% set the method to the default if not supplied
if nargin < 2
  method = 0;
end % if

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ECEF2LLA';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_ecef';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_ecef;
  
  estruct.variable(2).name = 'method';
  estruct.variable(2).req_dim = [1 1];
  estruct.variable(2).var = method;
  
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

if method == 0                 % use closed form approximation
  e2 = (RE^2 - b^2) / b^2;     % second numerical eccentricity (not e^2!)

  p = sqrt(x_ecef(:,1).^2 + x_ecef(:,2).^2);     % auxiliary quantity
  theta = atan2(RE .* x_ecef(:,3) , p .* b);     % auxiliary quantity

  % compute latitude (rad)
  lla(:,1) = atan2(x_ecef(:,3) + ...
                   e2 .* b .* sin(theta).^3, p - e2 .* RE .* cos(theta).^3);

  % compute longitude (rad)
  lla(:,2) = atan2(x_ecef(:,2), x_ecef(:,1));

  % compute the radius of curvature in the prime vertical
  n = RE^2 ./ sqrt(RE^2 * cos(lla(:,1)).^2 + b^2 * sin(lla(:,1)).^2);

  % compute height (m)
  lla(:,3) = p ./ cos(lla(:,1)) - n;   

else                      % use iterative method (method == 1)
  tol = 1e-14;            % set a tight tolerance on the change in lat (rad)
  max_iter = 15;          % set maximum iterations
  m = 1;                  % set first iteration
  dlat= 1;                % set the dlat value high for first iteration
  
  e2 = (RE^2 - b^2) / RE^2;                  % eccentricity^2
  p = sqrt(x_ecef(:,1).^2 + x_ecef(:,2).^2); % auxiliary quantity  
  lat = atan2(x_ecef(:,3), p .* (1 - e2));   % initial estimate of the latitude
  
  lla(:,2) = atan2(x_ecef(:,2), x_ecef(:,1)); % compute longitude
  
  % loop until the process has converged or maximum iterations are exceeded
  while (any(find(dlat > tol)) & m <= max_iter)
    % compute the radius of curvature in the prime vertical for this latitude
    n = RE^2 ./ sqrt(RE^2 * cos(lat).^2 + b^2 * sin(lat).^2); 
    
    lla(:,3) = (p ./ cos(lat)) - n;         % ellipsoidal height estimate 
    
    % new estiamte of lat
    newlat = atan2(x_ecef(:,3), p .* (1 - e2 .* (n ./ (n + lla(:,3)))));    
    
    dlat = abs(lat - newlat); % compute change in latitude for this iteration
    lat = newlat;             % update the latitude
    lla(:,1) = lat;           % set the return value to the last loop value
    m = m + 1;                % increase the loop counter
    if m > max_iter           % check against maximum iterations allowed
      fprintf('Warning: Maximum iterations (%d) exceeded in ECEF2LLA!\n', ...
               max_iter)
      lla(:,1) = lat;         % set the return value to the last loop value
      return
    end % if
  end % while

end % if method == 0 
  
%%%%% END ALGORITHM CODE %%%%%

% end ECEF2LLA      

