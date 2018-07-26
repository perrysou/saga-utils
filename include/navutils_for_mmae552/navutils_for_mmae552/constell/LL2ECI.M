function [x_eci] = ll2eci(x_ll, pos_eci, vel_eci)

% [x_eci] = ll2eci(x_ll, pos_eci, vel_eci)
%
% Function to convert vectors from Local Level (LL) to ECI coordinates.  
% This coordinate system is used as the local coordinate system for 
% satellites. The formulation must have a velocity vector to define the LL 
% y-axis. Inputs and outputs can be vectors. The LL frame used here is 
% defined as
%    z_hat_ll -> -R / |R|               (local down unit vector)
%    y_hat_ll -> -R cross V / |R x V|  (anti-orbit normal unit vector)
%    x_hat_ll -> y_hat cross z_hat
% 
% Input:
%   x_ll     - vector (LL) to be converted to ECI (m)
%               [xl, yl, zl], (nx3) 
%   pos_eci - satellite position (ECI) (m) 
%               [xe, ye, ze], (1x3 or nx3)
%   vel_eci - satellite velocity (ECI) (m/s)
%               [vxe, vye, vze], (1x3 or nx3)
% Output:
%   x_eci   - x_ll in eci frame [x, y, z] (nx3) (m)
%
% See also ECI2LL, LL2ECEF

% Written by: Maria J. Evans 9/3/97 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_eci=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(3,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on LL2ECI for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the 1x3 arguments to the full dimensions
if size(pos_eci,1) == 1 & size(pos_eci,2) == 3
  pos_eci = ones(size(x_ll,1),1) * pos_eci;
end
if size(vel_eci,1) == 1 & size(vel_eci,2) == 3
  vel_eci = ones(size(x_ll,1),1) * vel_eci;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'LL2ECI';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_ll';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_ll;
  
  estruct.variable(2).name = 'pos_eci';
  estruct.variable(2).req_dim = [1 3; 901 3];
  estruct.variable(2).var = pos_eci;
  
  estruct.variable(3).name = 'vel_eci';
  estruct.variable(3).req_dim = [1 3; 901 3];
  estruct.variable(3).var = vel_eci;
  
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

% Define the coordinate system

% Compute the magnitude of the radius vector
r_vect_mag = sqrt(pos_eci(:,1).^2 + pos_eci(:,2).^2 + pos_eci(:,3).^2);

% define the z-ll direction in ECI coordinates
z_hat_ll = -pos_eci ./ (r_vect_mag * ones(1,3));    

% find R x V and the magnitude of the R x V vector
r_cross_v = cross(pos_eci', vel_eci')';         
r_cross_v_mag = sqrt(r_cross_v(:,1).^2 + r_cross_v(:,2).^2 + r_cross_v(:,3).^2);

% define the y-ll direction in ECI coordinates
y_hat_ll = -r_cross_v ./ (r_cross_v_mag * ones(1,3));       

% define the x-ll direction in ECI coordinates
x_hat_ll = cross(y_hat_ll', z_hat_ll')';               

% Resize x_hat_ll, y_hat_ll, z_hat_ll from 1x3 to nx3, if necessary

if (size(x_hat_ll,1) == 1),     % x,y,z hat's are 1x3 --> convert to nx3
    x_hat_ll = ones(size(x_ll,1),1) * x_hat_ll;
    y_hat_ll = ones(size(x_ll,1),1) * y_hat_ll;
    z_hat_ll = ones(size(x_ll,1),1) * z_hat_ll;
end;

% write the ECI unit vectors in LL coordinates

x_hat_eci = [x_hat_ll(:,1),y_hat_ll(:,1),z_hat_ll(:,1)];
y_hat_eci = [x_hat_ll(:,2),y_hat_ll(:,2),z_hat_ll(:,2)];
z_hat_eci = [x_hat_ll(:,3),y_hat_ll(:,3),z_hat_ll(:,3)];

x_eci(:,1) = dot(x_ll',x_hat_eci')';    % x_eci component of x_ll
x_eci(:,2) = dot(x_ll',y_hat_eci')';    % y_eci component of x_ll
x_eci(:,3) = dot(x_ll',z_hat_eci')';    % z_eci component of x_ll
                                
%%%%% END ALGORITHM CODE %%%%%

% end of LL2ECI
