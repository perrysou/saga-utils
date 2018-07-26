function [x_body] = ned2body(x_ned, euler)

% [x_body] = ned2body(x_ned, euler);
%
% Function to transform a vector from the North-East-Down local-level frame to 
% a vehicle body frame.  The body frame is transformed from the NED frame by a
% 3-2-1 euler angle rotation sequence.  The rotation angles and axes of
% rotation are given by:
%
%     NED to NED' frame:    yaw about Down (3rd or z axis)
%     NED' to Body' frame:  pitch about E' (2nd or y axis)
%     Body' to Body frame:  roll about body-x axis (1st axis)
%
% Input:
%   x_ned - vector in NED frame (nx3) [North East Down]
%   euler - 3-2-1 Euler angle sequence from NED to body frame (radians) 
%            This input can be either (1x3) or (nx3), where the 
%            1st col. is yaw, 2nd col. is pitch, 3rd col. is roll.
%            If 1x3, then then the same Euler angles are used for each
%            vector. [yaw pitch roll] 
% Output:
%   x_body  - transformed vector [xb, yb, zb] (nx3)
%
% See also BODY2NED, ECEF2NED, NED2ECEF 

% Written by:  Maria J. Evans 2/13/97 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_body=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on NED2BODY for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Expland the euler variable to be the smae size as x_body if input is 1x3
if size(euler,1) == 1 & size(euler,2) == 3
  euler = ones(size(x_ned,1),1) * euler;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'NED2BODY';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_ned';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_ned;
  
  estruct.variable(2).name = 'euler';
  estruct.variable(2).req_dim = [901 3; 1 3];
  estruct.variable(2).var = euler;
  estruct.variable(2).type = 'ANGLE_RAD';
  
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

% allocate the output matrix
x_body = ones(size(x_ned,1),3) * inf;

% Compute the body vector a component at a time 

x_body(:,1) = cos(euler(:,1)).*cos(euler(:,2)).*x_ned(:,1) + ...
              sin(euler(:,1)).*cos(euler(:,2)).*x_ned(:,2) - ...
              sin(euler(:,2)).*x_ned(:,3);

x_body(:,2) = (-sin(euler(:,1)).*cos(euler(:,3)) + ...
                cos(euler(:,1)).*sin(euler(:,2)).*sin(euler(:,3))).*x_ned(:,1) + ...
               (cos(euler(:,1)).*cos(euler(:,3)) + ...
                sin(euler(:,1)).*sin(euler(:,2)).*sin(euler(:,3))).*x_ned(:,2) + ...
                cos(euler(:,2)).*sin(euler(:,3)).*x_ned(:,3);
             
x_body(:,3) = (sin(euler(:,1)).*sin(euler(:,3)) + ...
                cos(euler(:,1)).*sin(euler(:,2)).*cos(euler(:,3))).*x_ned(:,1) + ...
              (-cos(euler(:,1)).*sin(euler(:,3)) + ...
               sin(euler(:,1)).*sin(euler(:,2)).*cos(euler(:,3))).*x_ned(:,2) + ...
               cos(euler(:,2)).*cos(euler(:,3)).*x_ned(:,3);

%%%%% END ALGORITHM CODE %%%%%

% end NED2BODY
