function [x_ecef] = ned2ecef(x_ned, ref_lla)

% [x_ecef] = ned2ecef(x_ned, ref_lla);
%
% Function to convert vectors from NED (North-East-Down) coordinates to 
% ECEF (Earth-Centered-Earth-Fixed).
%
% Input:
%   x_ned   - NED position in meters [x_ned, y_ned, z_ned] (nx3)
%   ref_lla - reference lat, lon, and hgt for NED coordinate system base
%              lat and lon in rad, hgt in meters
%              ref_lla must have dimensions of 1x3 or nx3 [lat lon hgt]
%                 or 1x2 or nx2 [lat lon]
%              (Note:  hgt not used in this routine).
%              valid latitudes are -pi/2 -> pi/2
%              valid longitudes are -pi -> 2*pi
%
% Output:
%   x_ecef  - ECEF position in meters [x_ecef, y_ecef, z_ecef](nx3)
%
% See also ECEF2NED, ECEF2LLA, LLA2ECEF 

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
  fprintf('%s  See help on NED2ECEF for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Expland the ref_lla variable to be the smae size as x_ecef if input is 1x3
if size(ref_lla,1) == 1 & size(ref_lla,2) == 3
  ref_lla = ones(size(x_ned,1),1) * ref_lla;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'NED2ECEF';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_ned';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_ned;
  
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
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 
  
%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% allocate the output x_ecef vector and fill it with inf
% values for valid ref_lla will be filled in the following code
x_ecef = ones(size(x_ned)) * inf;

% compute the NED vector a component at a time             
x_ecef(:,1) = -x_ned(:,1) .* sin(ref_lla(:,1)) .* cos(ref_lla(:,2)) - ...
                     x_ned(:,2) .* sin(ref_lla(:,2)) - ...
                     x_ned(:,3) .* cos(ref_lla(:,1)) .* cos(ref_lla(:,2));
             
x_ecef(:,2) = -x_ned(:,1) .* sin(ref_lla(:,1)) .* sin(ref_lla(:,2)) + ...
                     x_ned(:,2)  .* cos(ref_lla(:,2)) - ...
                     x_ned(:,3) .* cos(ref_lla(:,1)) .* sin(ref_lla(:,2));
             
x_ecef(:,3) = x_ned(:,1) .* cos(ref_lla(:,1)) - x_ned(:,3) .* sin(ref_lla(:,1));

%%%%% END ALGORITHM CODE %%%%%

% end NED2ECEF
