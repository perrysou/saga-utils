function [x_ned] = ecef2ned(x_ecef, ref_lla)

% [x_ned] = ecef2ned(x_ecef, ref_lla);
%
% Function to convert vectors from ECEF (Earth-Centered-Earth-Fixed)
% to NED (North-East-Down) coordinates. 
%
% Input:
%   x_ecef  - ECEF position in meters (nx3) [x_ecef, y_ecef, z_ecef] 
%   ref_lla - reference lat, lon, and hgt for NED coordinate system base
%              lat and lon in rad.
%              ref_lla must have dimensions of (1x3 or nx3) or (1x2 or nx2)
%              [lat lon hgt]  or [lat lon]
%              (Note:  hgt is not used in this function)
%              valid latitudes are -pi/2 -> pi/2
%              valid longitudes are -pi -> 2*pi
%
% Output:
%   x_ned  - NED position in meters (nx3) [x_ned, y_ned, z_ned] 
%
% Note: This function is does not work at the poles.
%
% See also NED2ECEF, ECEF2LLA, LLA2ECEF 

% Written by: Jimmy Lamance 10/23/96 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

% Modified from 'Fundamentals of Astrodynamics',
% Bate, Mueller, and White page 101 conversion from IJK -> SEZ

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_ll=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECEF2NED for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Expland the ref_lla variable to be the smae size as x_ecef if input is 1x3
if (size(ref_lla,1) == 1 & size(ref_lla,2) == 3) | ...
   (size(ref_lla,1) == 1 & size(ref_lla,2) == 2) 
  ref_lla = ones(size(x_ecef,1),1) * ref_lla;
end

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.func_name = 'ECEF2NED';

estruct.variable(1).name = 'x_ecef';
estruct.variable(1).req_dim = [901 3];
estruct.variable(1).var = x_ecef;
  
estruct.variable(2).name = 'ref_lla';
estruct.variable(2).req_dim = [901 3; 1 3; 901 2; 1 2];
estruct.variable(2).var = ref_lla;
estruct.variable(2).type = 'LLA_RAD';
  
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% allocate the x_ned output matrix to all inf 
% valid x_ned's will be filled in the following code
x_ned = ones(size(x_ecef)) * inf;

% compute the NED vector a component at a time
x_ned(:,1) = -x_ecef(:,1) .* sin(ref_lla(:,1)) .* cos(ref_lla(:,2)) - ...
              x_ecef(:,2) .* sin(ref_lla(:,1)) .* sin(ref_lla(:,2)) + ...
              x_ecef(:,3) .* cos(ref_lla(:,1));
             
x_ned(:,2) = -x_ecef(:,1) .* sin(ref_lla(:,2)) + ...
              x_ecef(:,2)  .* cos(ref_lla(:,2));
             
x_ned(:,3) = -x_ecef(:,1) .* cos(ref_lla(:,1)) .* cos(ref_lla(:,2)) - ...
              x_ecef(:,2) .* cos(ref_lla(:,1)) .* sin(ref_lla(:,2)) - ...
              x_ecef(:,3) .* sin(ref_lla(:,1));             

%%%%% END ALGORITHM CODE %%%%%

% end ECEF2NED
