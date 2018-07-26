function I_pass = mask_vis(az, el, mask)

% I_pass = mask_vis(az, el, mask);
%
% Finds azimuth/elevation pairs passing the masking tests described 
% in the mask variable.
%
% Input:
%   az    - azimuth (rad) (nx1)
%            valid values are -2*pi -> 2*pi
%   el    - elevation (rad) (nx1)
%            valid values are -pi/2 -> pi/2
%   mask  - masking information (rad) (1x1, nx3, or nx4) (optional), default = 0;  
%        1x1 form is minimum elevation only [min_el]
%        nx3 form is min elevation and azimuth bounds [min_el start_az stop_az]
%        nx4 form is elevation and azimuth bounds [min_el max_el start_az stop_az]
%        Azimuth start and stop are assumed to be clockwise.
%        Examples:
%         minimum elevation mask of 5 degrees (rad) (1x1) 
%              mask = .0873   
%         minimum elevation and azimuth bound triples (nx3)
%              mask = [.0873   pi/2   pi;       % 5 deg min el, 90->180 azimuth
%                      .1745   pi    pi/4]      % 10 deg min el, 180->90 azimuth
%                                               % wraps through 0
%         elevation and azimuth bound quadruples
%              mask = [.0873  .5236  0    pi;    % 5->30 deg el, 0->180 azimuth
%                      .1745  pi/4   pi  2*pi]   % 10->45 deg el, 180->360 azimuth           
% Output:
%   I_pass - indices into the azimuth/elevation pairs that pass the masking test
%
% See also MASK_STA, NED2AZEL, ECEF2NED 

% Written by: Jimmy LaMance 11/8/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
lon=[]; lat=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on MASK_VIS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 3
  mask = 0;    % set mask to the default value of 0
end % if nargin < 3

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'MASK_VIS';

  % Develop the error checking structure with required dimension, matching
  % dinemsion flags, and input dinemsions.
  estruct.variable(1).name = 'az';
  estruct.variable(1).req_dim = [901 1];
  estruct.variable(1).var = az;
  estruct.variable(1).type = 'ANGLE_RAD';
  
  estruct.variable(2).name = 'el';
  estruct.variable(2).req_dim = [901 1];
  estruct.variable(2).var = el;
  estruct.variable(2).type = 'ELEVATION_RAD';
  
  estruct.variable(3).name = 'mask';
  estruct.variable(3).req_dim = [1 1; 902 3; 902 4];
  estruct.variable(3).var = mask;
  
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

% Verify that all azimuth definitions are from 0 -> 360 (not -180 -> 180)
I_az_negative = find(az < 0);                      % find all negative azimuths
az(I_az_negative) = az(I_az_negative) + 2 * pi;    % and convert to positive

if size(mask,2) == 1         % only an elevation mask is used
  I_pass = find(el >= mask); % find all the elevations above the elevation mask
  return                     % done! this is a really simple case

% else if we have an elevation and azimuth range to evaluate
elseif size(mask,2) == 3                   
  num_mask_pairs = size(mask,1);    % find out how many el/az pairs we have
  
  % define the return matrix of passing pairs as null ([])
  I_pass = [];                      
  
  for i = 1:num_mask_pairs 
  
    if mask(i,2) > mask(i,3)
      I_pass_new = find(el >= mask(i,1) & (az >= mask(i,2) | az <= mask(i,3)));    
    else  
      I_pass_new = find(el >= mask(i,1) & az >= mask(i,2) & az <= mask(i,3));
    end % if mask(i,2) > mask(i,3)
    
    I_pass = [I_pass; I_pass_new];
  end % for i = 1:num_mask_pairs  
  
  % sort back to the original ordering
  I_pass = sort(I_pass);

elseif size(mask,2) == 4                   
  num_mask_pairs = size(mask,1);    % find out how many el/az pairs we have
  
  % define the return matrix of passing pairs as null ([])
  I_pass = [];                      
  
  for i = 1:num_mask_pairs 
  
    if mask(i,3) > mask(i,4)
      I_pass_new = find(el >= mask(i,1) & el <= mask(i,2) & ...
                        (az >= mask(i,3) | az <= mask(i,4)));    
    else  
      I_pass_new = find(el >= mask(i,1) & el <= mask(i,2) & ...
                        az >= mask(i,3) | az <= mask(i,4));    
    end % if mask(i,3) > mask(i,4)
    
    I_pass = [I_pass; I_pass_new];
  end % for i = 1:num_mask_pairs  
  
  % sort back to the original ordering
  I_pass = sort(I_pass);

else
  fprintf('Unknown mask input format.  See help on mask_vis for details.\n');
  
end % if size(mask,2) == 1  

%%%%% END ALGORITHM CODE %%%%%

% end MASK_VIS
