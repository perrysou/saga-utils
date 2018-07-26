function [x_ll] = ecef2ll(x_vect, x_ecef, v_ecef)

% [x_ll] = ecef2ll(x_vect, x_ecef, v_ecef)
%
% Function to convert vectors from ECEF to Local Level (LL) coordinates.  This
% coordinate system is used as the local coordinate system for satellites.  
% The formulation must have a velocity vector to define the local-level y-axis. 
% Inputs and outputs can be vectors. The LL frame used here is defined as
%    z_hat_ll -> -R / |R|  (down)
%    y_hat_ll -> -R cross V / |R x V|  (anti-orbit-normal)
%    x_hat_ll -> y_hat cross z_hat
% 
% Input:
%   x_vect - vector (ECEF, meters) to be converted to LL in meters 
%             [x, y, z] (nx3)
%   x_ecef - satellite position (ECEF) in meters [xe, ye, ze]
%             1x3 or nx3
%   v_ecef - satellite velocity (ECEF) in m/s [vxe, vye, vze]
%             1x3 or nx3
%
% Output:
%   x_ll   - LL position in meters [xll, yll, zll] (nx3)
%
% Note: This function is does not work for a zero length velocity vector.
%
% See also LL2ECEF, ECEF2NED, NED2ECEF 

% Written by: Jimmy Lamance 11/7/96 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

% WGS-84 constants
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
x_ll=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(3,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ECEF2LL for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Expland the x_ecef and v_ecef variables to be the smae size as 
% x_body if they are input as 1x3
if size(x_ecef,1) == 1 & size(x_ecef,2) == 3
  x_ecef = ones(size(x_vect,1),1) * x_ecef;
end

if size(v_ecef,1) == 1 & size(v_ecef,2) == 3
  v_ecef = ones(size(x_vect,1),1) * v_ecef;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ECEF2LL';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'x_vect';
  estruct.variable(1).req_dim = [901 3];
  estruct.variable(1).var = x_vect;
  
  estruct.variable(2).name = 'x_ecef';
  estruct.variable(2).req_dim = [901 3; 1 3];
  estruct.variable(2).var = x_ecef;
  
  estruct.variable(3).name = 'v_ecef';
  estruct.variable(3).req_dim = [901 3; 1 3];
  estruct.variable(3).var = v_ecef;
  
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
r_vect_mag = sqrt(x_ecef(:,1).^2 + x_ecef(:,2).^2 + x_ecef(:,3).^2);

% define the z-ll direction in ECEF coordinates (the down vector)
z_hat_ll = -x_ecef ./ (r_vect_mag * ones(1,3));    

% add the Omega cross R term to the ECEF velocity so that we have
% inertial velocity expressed in the ECEF frame
e_rate_ecef = [0 0 EARTH_RATE]' * ones(1,size(x_ecef,1));
o_cross_r = cross(e_rate_ecef,x_ecef')';

% compute the ECI velocity in the ECEF frame
% (do a software no-no and use the same variable name for two different
%  physical quantities to save storage space)
v_ecef = v_ecef + o_cross_r;    % v_ecef is now really v_eci in ecef frame

% find R x V and the magnitude of the R x V vector in ECEF frame
r_cross_v = cross(x_ecef', v_ecef')';         
r_cross_v_mag = sqrt(r_cross_v(:,1).^2 + r_cross_v(:,2).^2 + r_cross_v(:,3).^2);

% define the y-ll direction in ECEF coordinates
y_hat_ll = -r_cross_v ./ (r_cross_v_mag * ones(1,3));        

% define the x-ll direction in ECEF coordinates
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

% end of ECEF2LL
