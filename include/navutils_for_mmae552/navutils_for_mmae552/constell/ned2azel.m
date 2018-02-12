function [az, el] = ned2azel(ned)

% [az, el] = ned2azel(ned);
%
% Function to convert a North, East, Down vector into azimuth and elevation.
%
%   Input:
%      ned - vector in north, east, and down coordinates (nx3)
%                (or could be used with any other coordinate system
%                 where north and east are two local-level unit vectors such
%                 that north x east = down, e.g., an aircraft navigation frame
%                 or a satellite local-level frame)
%   Output:
%      az  - azimuth (rad) (nx1):  rotation of vector in local North-East
%            plane, clockwise positive beginning at North,  0 <= az < 2*pi
%      el  - elevation (rad) (nx1):  angle of vector above North-East plane,
%            positive for a vector with a negative down component   (i.e., 
%                positive up component), -pi/2 <= el <= pi/2
%
% See also AZEL2NED, NED2BODY, BODY2NED, ECEF2NED, NED2ECEF 

% Written by: Jimmy LaMance 10/24/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
az=[]; el=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on NED2AZEL for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'NED2AZEL';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'ned';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = ned;
  
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

% azimuth is defined from north, therefore atan2(E/N)
az = atan2(ned(:,2),ned(:,1));                       % azimuth (rad)   

% Convert azimuths to be between 0 and 2*pi
I_neg_az = find(az < 0);
if ~isempty(I_neg_az)
  az(I_neg_az) = az(I_neg_az) + 2*pi;
end % if ~isempty(I_neg_az)

% elevation is defined from the local level (N-E) plane positive in the up 
% (-D) direction
north_east_mag = sqrt(ned(:,1).^2 + ned(:,2).^2);   % magnitude in the N-E plane
el = atan2(-ned(:,3), north_east_mag);              % elevation (rad)

%%%%% END ALGORITHM CODE %%%%%

% end NED2AZEL
