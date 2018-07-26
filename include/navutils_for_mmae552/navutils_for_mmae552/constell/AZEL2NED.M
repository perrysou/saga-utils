function [ned] = azel2ned(az, el)

% [ned] = azel2ned(az, el);
%
% Function to convert azimuth and elevation into an NED vector.
% Input:
%   az  - azimuth - rotation of vector about local vertical,
%         zero at North, pi/2 at East, etc. (rad) (nx1)
%         valid values are -2*pi -> 2*pi
%   el  - elevation: elevation angle from local-level,
%         pi/2 for straight up, -pi/2 for straight down (rad) (nx1)
%         valid values are -pi/2 -> pi/2
% Output:
%   ned - unit vector in north, east, and down coordinates (n x 3)
%
% See also NED2AZEL, BODY2NED

% Written by: Jimmy LaMance 4/2/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
ned=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on AZEL2NED for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'AZEL2NED';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'az';
  estruct.variable(1).req_dim = [901 1];
  estruct.variable(1).var = az; 
  estruct.variable(1).type = 'ANGLE_RAD';
  
  estruct.variable(2).name = 'el';
  estruct.variable(2).req_dim = [901 1];
  estruct.variable(2).var = el;
  estruct.variable(2).type = 'ELEVATION_RAD';
  
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

% initialize the output matrix
ned = ones(size(az,1),3) * inf;

cos_el = cos(el(:));      % compute just once here

ned(:,1) = cos_el.*cos(az(:));      % north component
ned(:,2) = cos_el.*sin(az(:));      % east component
ned(:,3) = -sin(el(:));             % down component

%%%%% END ALGORITHM CODE %%%%%

% end of AZEL2NED
