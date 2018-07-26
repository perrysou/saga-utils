function [trop_dry, trop_wet] = tropdlay(elev,trop_model)

% [trop_dry, trop_wet] = tropdlay(elev,trop_model);
%
% Computes the wet and dry troposphere delays using
% a Hopfield model.
%
% Input:
%   elev       - elevation angle to GPS satellites (rad) (nx1)
%                 Valid elevations are between 0 and pi/2.
%   trop_model - input for the troposphere model (optional)
%                 1x3 (or nx3) matrix in the form [p T e], where 
%                  p is the surface atmospheric pressure in mb,
%                  T is the surface temperature in degrees K, and
%                  e is the partial pressure of water vapor in mb.
%                 Default values are [1013.25 288.15 11.691].
% Output:
%   trop_dry   - Dry troposphere component (m) (nx1)
%   trop_wet   - Wet troposphere component (m) (nx1)
%
% See also PSEUDO_R 

% Written by: Jimmy LaMance 1/15/97
% Copyright (c) 1998 by Constell, Inc.

% Reference for Hopfield Model: GPS Theory and Practice, 
% Hofmann-Wellenhof, Lichtenegger, and Collins, pages 110-113. 
%
% Reference for default troposphere model parameters: 
% Collins, Langley, and LaMance, Limiting Factors in the Tropospheric
% Propagation Delay Error Modelling for GPS Airborne Navigation, 
% ION proceedings, June 1996, Boston, MA.

% functions called: ERR_CHK

% WGS-84 constants
RE = 6378137;     % WGS-84 value in meters

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
cphnew=[]; ambig=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on TROPDLAY for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Set the ionosphere model parameters to the default if not provided
if nargin < 2
  trop_model = [1013.25 288.15 11.691]; 
end % if nargin < 5

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'TROPDLAY';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'elev';
  estruct.variable(1).req_dim = [901 1];
  estruct.variable(1).var = elev;
  estruct.variable(1).type = 'ELEVATION_RAD_0';
  
  estruct.variable(2).name = 'trop_model';
  estruct.variable(2).req_dim = [1 3];
  estruct.variable(2).var = trop_model;
 
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% allocate the output data matrices and initialize to inf                                                                              
trop_dry = ones(size(elev,1),1) * inf;
trop_wet = ones(size(elev,1),1) * inf;

% convert the elevation angle to degrees
elev = elev * 180 / pi;

% define surface value constants    (Eqs 6.88 & 6.89)
c1 = 77.64;             % K / mb
c2 = -12.96;            % K / mb
c3 = 3.718e5;           % K^2 / mb

% compute dry and wet tropo refractive indices at 
% the surface (Eqs 6.88 & 6.89)
N_dry0 = c1 * trop_model(:,1) ./ trop_model(:,2);        
N_wet0 = c2 * trop_model(:,3) ./ trop_model(:,2) + ...
         c3 * trop_model(:,3) ./ trop_model(:,2).^2;

% define mean wet troposphere height (Eq 6.97)
h_mean_wet = ones(size(trop_model,1),1) * 11000;              % meters

% compute mean dry troposphere height  (Eq 6.91)
h_mean_dry = 40136 + 148.72 * (trop_model(:,2) - 273.16);     % meters

% compute the dry and wet tropospheric delays in meters (Eq 6.103)
dry_elev = sqrt(elev.^2 + 6.25) * pi / 180;   % convert it back to rad
wet_elev = sqrt(elev.^2 + 2.25) * pi / 180;   % convert it back to rad

trop_dry = (1.0e-6 / 5) .* N_dry0 .* h_mean_dry ./ sin(dry_elev);
trop_wet = (1.0e-6 / 5) .* N_wet0 .* h_mean_wet ./ sin(wet_elev);

%%%%% END ALGORITHM CODE %%%%%

% end of TROPDLAY
  
