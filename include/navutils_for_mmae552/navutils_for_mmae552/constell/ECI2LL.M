function [x_ll] = eci2ll(x_vect, x_eci, v_eci)

% [x_ll] = eci2ll(x_vect, x_eci, v_eci);
%
% Function to convert vectors from ECI to Local Level (LL) coordinates.  This
% coordinate system is used as the local coordinate system for satellites.  
% The formulation must have a velocity vector to define the local-level y-axis. 
% Inputs and outputs can be vectors. The LL frame used here is defined as
%    z_hat_ll -> -R / |R|  (down)
%    y_hat_ll -> -R cross V / |R x V|  (anti-orbit-normal)
%    x_hat_ll -> y_hat cross z_hat
% 
% Input:
%   x_vect - vector (ECI, meters) to be converted to LL in meters 
%             [x, y, z] (nx3)
%   x_eci - satellite position (ECI) in meters [xi, yi, zi]
%             1x3 or nx3
%   v_eci - satellite velocity (ECI) in m/s [vxi, vyi, vzi]
%             1x3 or nx3
%
% Output:
%   x_ll   - LL position in meters [xll, yll, zll] (nx3)
%
% Note: This function is does not work for a zero length velocity vector.
%
% See also ECEF2LL, LL2ECI, LL2ECEF

% Written by: Maria J. Evans 9/3/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_ll=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(3,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECI2LL for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Expland the euler variable to be the smae size as x_body if input is 1x3
if size(x_eci,1) == 1 & size(x_eci,2) == 3
  x_eci = ones(size(x_vect,1),1) * x_eci;
end

if size(v_eci,1) == 1 & size(v_eci,2) == 3
  v_eci = ones(size(x_vect,1),1) * v_eci;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ECI2LL';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_vect';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_vect;
  
  estruct.variable(2).name = 'x_eci';
  estruct.variable(2).req_dim = [901 3; 1 3];
  estruct.variable(2).var = x_eci;
  
  estruct.variable(3).name = 'v_eci';
  estruct.variable(3).req_dim = [901 3; 1 3];
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

% Define the coordinate system

% Compute the magnitude of the radius vector
r_vect_mag = sqrt(x_eci(:,1).^2 + x_eci(:,2).^2 + x_eci(:,3).^2);

% define the z-ll direction in ECI coordinates (the down vector)
z_hat_ll = -x_eci ./ (r_vect_mag * ones(1,3));    

% find R x V and the magnitude of the R x V vector in ECI frame
r_cross_v = cross(x_eci', v_eci')';         
r_cross_v_mag = sqrt(r_cross_v(:,1).^2 + r_cross_v(:,2).^2 + r_cross_v(:,3).^2);

% define the y-ll direction in ECI coordinates
y_hat_ll = -r_cross_v ./ (r_cross_v_mag * ones(1,3));        

% define the x-ll direction in ECI coordinates
x_hat_ll = cross(y_hat_ll', z_hat_ll')';               

% Resize x_hat_ll, y_hat_ll, z_hat_ll from 1x3 to nx3, if necessary

if (size(x_hat_ll,1) == 1),     % x,y,z hat's are 1x3 --> convert to nx3
    x_hat_ll = ones(size(x_vect,1),1) * x_hat_ll;
    y_hat_ll = ones(size(x_vect,1),1) * y_hat_ll;
    z_hat_ll = ones(size(x_vect,1),1) * z_hat_ll;
end;

x_ll(:,1) = dot(x_vect', x_hat_ll')'; % compute local_level x-axis output
x_ll(:,2) = dot(x_vect', y_hat_ll')'; % compute local_level y-axis output
x_ll(:,3) = dot(x_vect', z_hat_ll')'; % compute local_level z-axis output

%%%%% END ALGORITHM CODE %%%%%

% end of ECI2LL
