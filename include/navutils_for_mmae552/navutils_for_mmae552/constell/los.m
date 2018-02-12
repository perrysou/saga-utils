function [t_los,los_vect,los_indices,obscure_info] = los(t1,x1,t2,x2, ...
                                                      earth_model_flag);

% [t_los,los_vect,los_indices,obscure_info] = los(t1,x1,t2,x2,earth_model_flag)
% 
% Compute line-of-sight (LOS) vectors from a group of objects to another
% group of objects, each identified by object number within the group.  For 
% example, this function allows computation of LOS vectors from a satellite
% (group #1) to the GPS constellation of satellites (group #2).  Group #2
% could also be a group of ground stations.  This function is used for LOS
% computation for all ground assets, satellites, and constellation combinations.
%
% Input:
%   t1  - GPS time vector for group #1 objects (mx2) [GPS_week GPS_sec]
%          valid GPS_week values are 1-3640 (years 1980-2050)
%          valid GPS_sec values are 0-604799
%   x1  - positions for group #1 vehicles (mx3) [x y z] or (mx4) [veh_num1 x y z]
%   t2  - GPS time vector for group #2 objects (nx2) [GPS_week GPS_sec]
%          valid GPS_week values are 1-3640 (years 1980-2050)
%          valid GPS_sec values are 0-604799
%   x2  - positions for group #2 vehicles (nx3) [x y z] or (nx4) [veh_num2 x y z], 
%          veh_num is a vehicle identification number, not required for either
%          group #1 or group #2.
%   earth_model_flag - 0 for spherical earth, 1 for oblate earth (1x1).
%                      Default = spherical earth. (Optional). Used only if
%                      obscure_info is requested as output. The tangent of the los_vect
%                      to the earth is computed. Note: For oblate earth, input vectors
%                      x1 and x2 must be in ECEF coordinates.                    
%                      Not used for visibility from a ground-based site.
%
%   Note: If a single time for a given object number is given for either input set
%         (e.g. x1 or x2 data), it is assumed to be an object fixed to the surface
%         of the Earth.  This objects position will be fixed and LOS computed
%         for each of the times in the other object (e.g. if the ground station
%         is entered in x1 data, then the times from x2 will be used for the
%         LOS computations).  The Earth fixed objects are assumed to be in ECEF
%         coordinates.  Care should be taken when inputing Earth fixed objects
%         in ECI or other non-ECEF coordinate frames.
%   Note: Group #1 and group #2 positions must be in the same reference system 
%          (e.g. ECEF or ECI).
% Output:
%   t_los        - GPS time vector corresponding to the LOS (kx2) [GPS_week GPS_sec]
%   los_vect          - line of sight vector at t_los from group #1 object to 
%                  group #2 object (kx3).  LOS are computed only 
%                  for matching times in t1 and t2.
%   los_indices  - indices to obtain relationship between output LOS vectors and
%                  the input positions [x1_ind, x2_ind] (kx2)
%   obscure_info - contains information needed to determine whether the earth
%                  obscures the line-of-sight (kx3) [tangent_altitude, alpha, beta].
%                  Alpha is the angle between the x1 and x2 vectors. 
%                  Beta is the angle between the x1 vectors and the radius to the
%                  tangent point of the los vectors.
%                  An observation is obscured if the tangent vector magnitude is
%                  below the users tolerance, and alpha is greater than beta. 
%
% Note: The los_indices can be used to track which velocities or other
%       data associated with the vehicle trajectory go with each line of sight.
%       For example, to rotate the ECEF LOS vector to a satellite local
%       level coordinate system via ECEF2LL, the position and velocity of the
%       vehicle are required.  The ECEF2LL function call would be the following
%       assuming the original vehicle position and velocity matrices were x1 & v1
%         [x_ll] = ecef2ll(los_vect, x1(los_indices(:,1),:), v1(los_indices(:,1),:));
%       See help on ECEF2LL for details about these inputs.
%
% See also PROPGEPH

% Written by: Jimmy LaMance 10/24/96
% Modified by: Maria Evans, May 98
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPST2SEC, REMAP, SEC2GPST, TANELLIP

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_los = []; los_vect = []; veh_num1 = []; veh_num2 = []; x1_out = []; x2_out = [];
x1_ind = []; x2_ind = [];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(4,5,nargin);
if ~isempty(msg)
  fprintf('%s  See help on LOS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

estruct.func_name = 'LOS';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 't1';
estruct.variable(1).req_dim = [901 2];
estruct.variable(1).var = t1;
estruct.variable(1).type = 'GPS_TIME';
  
estruct.variable(2).name = 'x1';
estruct.variable(2).req_dim = [901 3; 901 4];
estruct.variable(2).var = x1;
  
estruct.variable(3).name = 't2';
struct.variable(3).req_dim = [902 2];
estruct.variable(3).var = t2;
estruct.variable(3).type = 'GPS_TIME';
  
estruct.variable(4).name = 'x2';
estruct.variable(4).req_dim = [902 3; 902 4];
estruct.variable(4).var = x2;
  
if nargin==5,
  estruct.variable(4).name = 'earth_model_flag';
  estruct.variable(4).req_dim = [1 1];
  estruct.variable(4).var = earth_model_flag;
else,
  earth_model_flag = 0;
end;
 
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
RE = 6378137.0;

% If x1 is nx3, resize to be nx4 with 1's in the satellite number column
if size(x1,2) == 3 
  x1(:,2:4) = x1;
  x1(:,1) = ones(size(x1(:,1)));
end % size(x1,2) == 3

if size(x2,2) == 3 
  x2(:,2:4) = x2;
  x2(:,1) = ones(size(x2(:,1)));
end % size(x2,2) == 3

% Convert to linear time in seconds past the GPS epoch                                      
t1 = gpst2sec(t1);
t2 = gpst2sec(t2);

% Find all of the single times/position for the group #1 and group #2 input data.
% These are assumed to be fixed to the Earth and the ECEF positions (assumed ECEF)
% will be duplicated for all of the times in the other group set.
% Check the size of the t1 input, if it is 1x1, then expand it to be
% the same times as t2 and assume the position is constant
% Find the unique input times for the combined x1 and x2 data
t_all = unique([t1; t2]);

% Find the unique vehicle numbers in data set 1
vn1 = unique(x1(:,1));

% Set up indices indicating if this is a ground station in the group.  Used later
% in re-establishing the output indices.
ind_x1 = zeros(length(vn1),1);
sta_num_x1 = zeros(length(vn1),1);

x1_org = x1;
x2_org = x2;

% Loop over the vehicle #1 data to find the unique times
for i = 1:length(vn1)
  I = find(x1_org(:,1) == vn1(i));
  
  % Find the unique times for this data
  tu = unique(t1(I));
  
  % If this is a single time, expand it to be the size of t_all
  if length(tu) == 1
    % Start by removing the single time point.  Prevents repetition of times
    I_t1 = find(x1(:,1) ~= vn1(i));
    t1 = t1(I_t1);
    x_to_duplicate = x1_org(I(1),:);
    x1 = x1(I_t1,:);
    
    % Add all of the t_all times
    t1 = [t1; t_all];
    x1 = [x1; ones(length(t_all),1) * x_to_duplicate];
    ind_x1(i) = I(1);  
    sta_num_x1(i) = x1_org(I(1),1);

  end % if length(tu) == 1
  
%  clear I tu x_to_duplicate
end % for i = 1:length(vn1)

% Repeat the process for the group #2 data
% Find the unique vehicle numbers in data set 1
vn2 = unique(x2(:,1));
ind_x2 = zeros(length(vn2),1);
sta_num_x2 = zeros(length(vn2),1);

% Loop over the vehicle #1 data to find the unique times
for i = 1:length(vn2)
  I = find(x2_org(:,1) == vn2(i));
  
  % Find the unique times for this data
  tu = unique(t2(I));
  
  % If this is a single time, expand it to be the size of t_all
  if length(tu) == 1
    % Start by removing the single time point.  Prevents repetition of times
    I_t2 = find(x2(:,1) ~= vn2(i));
    t2 = t2(I_t2);
    x_to_duplicate = x2_org(I(1),:);
    x2 = x2(I_t2,:);
    
    % Add all of the t_all times
    t2 = [t2; t_all];
    x2 = [x2; ones(length(t_all),1) * x_to_duplicate];
    ind_x2(i) = I(1);
    sta_num_x2(i) = x2_org(I(1),1);
  end % if length(tu) == 1
  
  clear I tu x_to_duplicate
end % for i = 1:length(vn2)

% Check to see if the data need to be time sorted  
dt1 = diff(t1);
dt2 = diff(t2);

if any(dt1 < 0);
  [t1, I_sort] = sort(t1);
  x1 = x1(I_sort,:);
end
if any(dt2 < 0);
  [t2, I_sort] = sort(t2);
  x2 = x2(I_sort,:);
end

% Find the intersection of the 2 time sets.  These are the unique
% intersection points with duplicates deleted.  The indices point to 
% the last of a run of common points (eg t1 = [1 1 1] -> I_t1 = 3)
[t_intersect, I_t1, I_t2] = intersect(t1,t2);

% Make sure that there are some common time intersections
if isempty(t_intersect)
  fprintf('No common times in t1 and t2 data in LOS.\n')
  return
end

% Set up indices to all the common times.  This will also include
% the non-common times at this point.  The non-common times will
% be removed later.
t1_all_index = zeros(size(t1));
t2_all_index = zeros(size(t2));

num_t1_intersect = length(I_t1);
num_t2_intersect = length(I_t2);

% Remove the last intersection point from both index data sets
I_t1 = I_t1(1:end-1);
I_t2 = I_t2(1:end-1);

% Insert a 1 at all of the intersection time changes.
t1_all_index(I_t1+1) = ones(size(I_t1));
t2_all_index(I_t2+1) = ones(size(I_t2));

% Put a 1 at the starting index
t1_all_index(1) = 1;
t2_all_index(1) = 1;           

% Use the cumsum function to create the indices to the intersection of the 
% times.  This will results in matrices like [1 1 1 1 2 2 2 2... n n n].
t1_all_index = cumsum(t1_all_index);
t2_all_index = cumsum(t2_all_index); 

% Set up another set of time matrices that are the same size as the initial
% matrices, but only contain times that are the intersection of the two time
% matrices.  This will be used to eliminate the non-intersection times.
t1_new = t_intersect(t1_all_index);
t2_new = t_intersect(t2_all_index);

% Find all of the valid
I_t1_keep = find(t1 == t1_new);
I_t2_keep = find(t2 == t2_new);

% Eliminate time tags for data not used in the comparison
t1 = t1(I_t1_keep);
x1 = x1(I_t1_keep,:);
t2 = t2(I_t2_keep);
x2 = x2(I_t2_keep,:);
t1_index = t1_all_index(I_t1_keep);
t2_index = t2_all_index(I_t2_keep);

% Find all of the unique satellite/station numbers for both data sets.  There will
% be a loop over the smaller of the two sets.
prn_unique_1 = unique(x1(:,1));
t_unique_1 = unique(t1);
prn_unique_2 = unique(x2(:,1));
t_unique_2 = unique(t2);

num_prn1 = length(prn_unique_1);
num_t1 = length(t_unique_1);
num_prn2 = length(prn_unique_2);
num_t2 = length(t_unique_2);

% Switch group #1 and group #2 if required to make group #1 have the fewest 
% number of vehicles.
swap_12 = 0;               % flag indicating if the swap has occured
if num_prn1 > num_prn2
  t2_temp = t2;
  x2_temp = x2;
  t2_index_temp = t2_index;
  t2 = t1;
  x2 = x1;
  t2_index = t1_index;
  t1 = t2_temp;
  x1 = x2_temp;       
  t1_index = t2_index_temp;
  
  prn_unique_1 = unique(x1(:,1));
  t_unique_1 = unique(t1);
  prn_unique_2 = unique(x2(:,1));
  t_unique_2 = unique(t2);

  num_prn1 = length(prn_unique_1);
  num_t1 = length(t_unique_1);
  num_prn2 = length(prn_unique_2);
  num_t2 = length(t_unique_2);  
  
  I_t1_keep_temp = I_t1_keep;
  I_t1_keep = I_t2_keep;
  I_t2_keep = I_t1_keep_temp;

  swap_12 = 1;                       % set the swap flag
end % if num_prn1 > num_prn2

% Remap the PRN 2 data to be indexed as 1, 2, ... num_prn2
[prn2_index] = remap(x2(:,1));

% Create a matrix filled with inf the size of the satellite 2 data
% assuming that each colunm is a time and the rows are satellite numbers
d2_prn = ones(num_prn2,num_t2) * inf;      % vehicle number
d2_x = ones(num_prn2,num_t2) * inf;        % X-component
d2_y = ones(num_prn2,num_t2) * inf;        % Y-component
d2_z = ones(num_prn2,num_t2) * inf;        % Z-component
d2_t = ones(num_prn2,num_t2) * inf;        % time
  
% Fill each of the times and available PRN data with ones for satellite 2.
% This index is used to fill in the data for satellite/constellation #2.
I_d2 = sub2ind(size(d2_x),prn2_index',t2_index);

% Fill in the X, Y, and Z data for the group #2 assets
d2_prn(I_d2) = x2(:,1);
d2_x(I_d2) = x2(:,2);
d2_y(I_d2) = x2(:,3);
d2_z(I_d2) = x2(:,4);
d2_t(I_d2) = t2;

% Check that the matrix dimension for the x2 data match.  This can cause
% a problem when there is repeated time data points without the appropriate
% satellite numbers.
if size(I_t2_keep,1) ~= size(d2_x,2) *  size(d2_x,1)
  fprintf('Problem detected in LOS.\n')
  fprintf('Verify that the any constellation position data has the satellite\n')
  fprintf('numbers in addition to the position information.  Also verify that\n')
  fprintf('there are no repeated time/position/satellite number combinations.\n')
  fprintf('Returning from LOS with empty output matrices.\n\n');
  dbstack
  fprintf('\n\n')                                            
  return
end % if size(I_t2_keep,1) ~= size(d2_x,2)

% Initialize the output variable 
los_vect = [];
t_los_sec = [];
veh_num1 = [];
veh_num2 = [];
x1_out = [];
x2_out = [];
x1_ind = [];
x2_ind = [];

% Loop over the vehicle numbers in group #1
for i = 1:num_prn1
  % Find all of the time #1 where this PRN is used.  This will have a single
  % index for each time.
  this_prn = prn_unique_1(i);
  I_times_1 = find(x1(:,1) == this_prn);
                                       
  % Get the indices to the time when this PRN is available
  d1_index = t1_index(I_times_1)';                         
  
  % Create a vector with infinity the size of the whole common time series
  d1_row_one_x = ones(1,size(d2_x,2)) * inf;
  d1_row_one_y = ones(1,size(d2_x,2)) * inf;
  d1_row_one_z = ones(1,size(d2_x,2)) * inf;
  
  % Fill in x, y, and z position data at time when this PRN is used
  d1_row_one_prn(d1_index) = x1(I_times_1,1);
  d1_row_one_x(d1_index) = x1(I_times_1,2);
  d1_row_one_y(d1_index) = x1(I_times_1,3);
  d1_row_one_z(d1_index) = x1(I_times_1,4);
  
  % Make this the same size as the d1 (group #2) data
  d1_prn = ones(size(d2_x,1),1) * d1_row_one_prn;
  d1_x(i,:,:) = ones(size(d2_x,1),1) * d1_row_one_x;
  d1_y(i,:,:) = ones(size(d2_x,1),1) * d1_row_one_y;
  d1_z(i,:,:) = ones(size(d2_x,1),1) * d1_row_one_z;   

  d2_x_all(i,:,:) = d2_x;
  d2_y_all(i,:,:) = d2_y;
  d2_z_all(i,:,:) = d2_z;   
  t_los_sec(i,:,:) = d2_t;

  veh_num1(i,:,:) = d1_prn;
  veh_num2(i,:,:) = d2_prn;
  x1_ind(i,:,:) = ones(size(d2_x,1),1) * I_t1_keep(I_times_1)'; 
  x2_ind(i,:,:) = reshape(I_t2_keep,size(d2_x)); 
  
end

% Compute the LOS for this group #1 vehicle
los_x = d2_x_all - d1_x;
los_y = d2_y_all - d1_y;
los_z = d2_z_all - d1_z;

% Remove the LOS that are INF (no group #2 or group #1 data)
I_valid_los = find(isinf(los_x) == 0);

num_los = length(I_valid_los);    

% Reshape the output to be vectors while eliminating the invalid points
los_x = reshape(los_x(I_valid_los),num_los,1);
los_y = reshape(los_y(I_valid_los),num_los,1);
los_z = reshape(los_z(I_valid_los),num_los,1);
t_los_sec = reshape(t_los_sec(I_valid_los),num_los,1);
veh_num1 = reshape(veh_num1(I_valid_los),num_los,1);  
veh_num2 = reshape(veh_num2(I_valid_los),num_los,1);  
  
x1_ind = reshape(x1_ind(I_valid_los),num_los,1);
x2_ind = reshape(x2_ind(I_valid_los),num_los,1);

d1_x = reshape(d1_x(I_valid_los),num_los,1);
d1_y = reshape(d1_y(I_valid_los),num_los,1);
d1_z = reshape(d1_z(I_valid_los),num_los,1);

d2_x_all = reshape(d2_x_all(I_valid_los),num_los,1);
d2_y_all = reshape(d2_y_all(I_valid_los),num_los,1);
d2_z_all = reshape(d2_z_all(I_valid_los),num_los,1);

los_vect = [los_x los_y los_z];
x1_out = [d1_x d1_y d1_z];
x2_out = [d2_x_all d2_y_all d2_z_all];
  
% See if the group #1 and group #2 data were swapped, if so swap them back
if swap_12 == 1                                                           
  los_vect = -los_vect;
  veh_num1_temp = veh_num1;
  veh_num1 = veh_num2;
  veh_num2 = veh_num1_temp;
  x1_out_temp = x1_out;
  x1_out = x2_out;
  x2_out = x1_out_temp;
  x1_ind_temp = x1_ind;
  x1_ind = x2_ind;
  x2_ind = x1_ind_temp;

end % if swap_12 == 1

% Clean up x1 and x2 indices if single ground stations time/positions were input
I_x1 = find(ind_x1 ~= 0);
if any(I_x1)
  for i = 1:length(I_x1)
    this_sta_num = sta_num_x1(I_x1(i));
  
    I_replace = find(this_sta_num == veh_num1);
    
    % Replace the corresponding indices in the output 
    x1_ind(I_replace) = ones(size(I_replace)) * ind_x1(i);
  end % for i = 1:length(I_x1)
end % if any(I_x1)

% Clean up x1 and x2 indices if single ground stations time/positions were input
I_x2 = find(ind_x2 ~= 0);
if any(I_x2)
  
  for i = 1:length(I_x2)
    this_sta_num = sta_num_x2(I_x2(i));
  
    I_replace = find(this_sta_num == veh_num2);
    
    % Replace the corresponding indices in the output 
    x2_ind(I_replace) = ones(size(I_replace)) * ind_x2(i);
  end % for i = 1:length(I_x2)
end % if any(I_x2)

% Get time back into full GPS time
t_los = sec2gpst(t_los_sec);

% Fill indices
los_indices = [x1_ind, x2_ind];

if nargout==4,  % Compute the Earth obscuring information if requested
  % Call tanellip to compute obscuring info
  obscure_info = tanellip(earth_model_flag,x1_out,x2_out,los_vect);
end;  % if nargout==4,

%%%%% END ALGORITHM CODE %%%%%

% end LOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  TANELLIP function within LOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obscure_info] = tanellip(earth_model_flag,x1_out,x2_out,los_vect);
% Compute tangent to ellipsoid
RE = 6378137.0;
seminor = 6356752.314;         % ellipsiod semi-minor axis


%%%%% BEGIN ALGORITHM CODE %%%%%
if earth_model_flag==0,  % Use spherical earth
  % Compute the radius to the tangent point to the Earth's surface
  % First compute the angle between the los vector and x1_out
  [x1_norm, x1_mag] = normvect(x1_out);
  [los_norm, los_mag] = normvect(los_vect);
  [x2_norm, x2_mag] = normvect(x2_out);
  
  % Set the initial values for the tangent altitudes to the target vehicles
  tangent_altitude = x2_mag - RE;
  
  delta = pi - acos(dot(x1_out',los_vect')'./(x1_mag.*los_mag));
  beta = zeros(size(delta));
  below_90 = find(delta < pi/2);
  
  if any(below_90),
    beta(below_90) = (pi/2 - delta(below_90));
    tangent_altitude(below_90) = cos(beta(below_90)).*x1_mag(below_90) - RE;
  end;
  % Cases with delta above pi/2 are clearly not obscured by the earth,
  % so set beta to pi so that alpha is never greater than beta
  above_90 = find(delta > pi/2);
  if any(above_90),
    beta(above_90) = ones(length(above_90),1)*pi;
  end;
  equal_90 = find(delta == pi/2);
  if any(equal_90),
    beta(equal_90) = zeros(length(equal_90),1);
  end;

  alpha = acos(dot(x1_out',x2_out')'./(x1_mag.*x2_mag));
  obscure_info = [tangent_altitude, alpha, beta];

else,  % Oblate Earth model
  % Iterate to compute the lat and lon of the tangent altitude.
  % Use FMIN function to compute the minimum

  % First compute the angle between the los vector and x1_out
  [x1_norm, x1_mag] = normvect(x1_out);
  [los_norm, los_mag] = normvect(los_vect);
  [x2_norm, x2_mag] = normvect(x2_out);

  % Set the initial tangent altitude to be the altitude of the x1_out positions.
  % The tangent altitude that correspond to delta angles less than 90
  % will be recomputed below.
  x1_lla = ecef2lla(x1_out,1);
  tangent_altitude = x1_lla(:,3);

  % Compute the angle, delta, between the x1 position and the LOS vector
  delta = pi - acos(dot(x1_out',los_vect')'./(x1_mag.*los_mag));
  
  % Find the delta angles less than 90 degrees
%  below_90 = find(delta < pi/2);

%  if any(below_90),
    % these are potentially cases which my be obscured by the oblate earth
    % use FMIN to find the minimum altitude above the ellipsoid
%    [los_tangent_length] = fminv('ellipalt',zeros(size(los_mag(below_90))),...
%                               los_mag(below_90),[0,1.e-12], ...
%                               x1_out(below_90,1), x1_out(below_90,2),...
%                               x1_out(below_90,3), los_norm(below_90,1),...
%                               los_norm(below_90,2),los_norm(below_90,3));
%    rt(:,1) = x1_out(below_90,1) + los_norm(below_90,1).*los_tangent_length;
%    rt(:,2) = x1_out(below_90,2) + los_norm(below_90,2).*los_tangent_length;
%    rt(:,3) = x1_out(below_90,3) + los_norm(below_90,3).*los_tangent_length;
    
%    [rt_norm, rt_mag] = normvect(rt);
%    [lla] = ecef2lla(rt);
    
%    tangent_altitude(below_90) = lla(:,3);
%    beta(below_90,1) = asin(sin(delta(below_90)).*los_tangent_length./rt_mag);
%  end;
  
  % Cases with delta above pi/2 are clearly not obscured by the earth,
  % so set beta to pi so that alpha is never greater than beta
%  above_90 = find(delta > pi/2);
%  if any(above_90),
%    beta(above_90,1) = ones(length(above_90),1)*pi;
%  end;
%  equal_90 = find(delta == pi/2);
%  if any(equal_90),
%    beta(equal_90,1) = zeros(length(equal_90),1);
%  end;

  [los_tangent_length] = fminv('ellipalt',zeros(size(los_mag)),...
                             los_mag,[0,1.e-12], ...
                             x1_out(:,1), x1_out(:,2),...
                             x1_out(:,3), los_norm(:,1),...
                             los_norm(:,2),los_norm(:,3));
  rt(:,1) = x1_out(:,1) + los_norm(:,1) .* los_tangent_length;
  rt(:,2) = x1_out(:,2) + los_norm(:,2) .* los_tangent_length;
  rt(:,3) = x1_out(:,3) + los_norm(:,3) .* los_tangent_length;
    
  [rt_norm, rt_mag] = normvect(rt);
  [lla] = ecef2lla(rt);
    
  tangent_altitude(:) = lla(:,3);
  beta = asin(sin(delta).*los_tangent_length ./ rt_mag);

  alpha = acos(dot(x1_out',x2_out')'./(x1_mag.*x2_mag));

  obscure_info = [tangent_altitude, alpha, beta];
end;

%%%%% END ALGORITHM CODE %%%%%

% end TANELLIP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  REMAP function within LOS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_remap = remap(d);

% d_remap = remap(d);
%
% Remaps a vector of data with repeated indices to the values 1 through
% the number of unique entires in d.  This function can be used to take
% a set of satellite numbers and reindex them to be satellites 1 through the
% total number of unique satellites without skipping numbers.  This is a 
% support function for COMP_LOS.
%
% Input: d - vector of numbers to be indexed as as 1, 2, ...n where n is
%            the number of unique inputs.
% 
% Output: d_remap - indexed version of the input, d, where the unique numbers
%                   have been replace with 1 through the maximum number of 
%                   unique numbers

% Written by: Jimmy LaMance 4/14/98
% Copyright (c) 1998 by Constell, Inc.

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
d_remap=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on REMAP for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

estruct.func_name = 'REMAP';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'd';
estruct.variable(1).req_dim = [901 1];
estruct.variable(1).var = d;

% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% Start by sorting the input data
[d_sort, ndx_sort] = sort(d);

% Find all of the unique data points
[d_unique, ndx_unique] = unique(d_sort);

% Remove the last unique index
ndx_unique = ndx_unique(1:end-1);

% Build a matrix with zeros the length of d and put 1's at the changes
z = zeros(size(d));
z(ndx_unique+1) = ones(size(ndx_unique));

% Add a one to start the z vector
z(1) = 1;

% Use the cumsum function to get to 1, 2, 3, ...
z_cum = cumsum(z);     

% Resort the remapped data to the original indices
d_remap(ndx_sort) = z_cum;

%%%%% End of REMAP %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  FMINV function within LOS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xf,options] = fmin(funfcn,ax,bx,options,varargin)
%FMIN   Minimize function of one variable.
%   X = FMIN('F',x1,x2) attempts to return a value of x which is a local 
%   minimizer of F(x) in the interval x1 < x < x2.  'F' is a string 
%   containing the name of the objective function to be minimized.
%
%   X = FMIN('F',x1,x2,OPTIONS) uses a vector of control parameters.
%   If OPTIONS(1) is positive, intermediate steps in the solution are
%   displayed; the default is OPTIONS(1) = 0.  OPTIONS(2) is the termination
%   tolerance for x; the default is 1.e-4.  OPTIONS(14) is the maximum
%   number of function evaluations; the default is OPTIONS(14) = 500.  
%   The other components of OPTIONS are not used as input control 
%   parameters by FMIN.  For more information, see FOPTIONS.
%
%   X = FMIN('F',x1,x2,OPTIONS,P1,P2,...) provides for additional
%   arguments which are passed to the objective function, F(X,P1,P2,...)
%
%   [X,OPTIONS] = FMIN(...) returns a count of the number of steps
%   taken in OPTIONS(10).
%
%   Examples
%       fmin('cos',3,4) computes pi to a few decimal places.
%       fmin('cos',3,4,[1,1.e-12]) displays the steps taken
%       to compute pi to about 12 decimal places.
%
%   See also FMIN, FMINS.

%   Reference: "Computer Methods for Mathematical Computations",
%   Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.

% Taken from original MATLAB function FMIN.
% Vectorized by Jimmy LaMance, June 1998
% Constell, Inc. 

% initialization
if nargin<4, options = []; end
options = foptions(options);
tol = options(2);
if (~options(14))
    options(14) = 500; 
  end
if nargin < 3
  error('fmin requires three input arguments.')
end

header = ' Func evals     x           f(x)         Procedure';
step='       initial';

% Convert to inline function as needed.
funfcn = fcnchk(funfcn,length(varargin));

num = 1;
seps = sqrt(eps);
c = 0.5*(3.0 - sqrt(5.0));
a = ax; b = bx;
v = a + c*(b-a);
w = v; xf = v;
d = 0.0; e = 0.0;
x= xf; fx = feval(funfcn,x,varargin{:});  
% fmin_data = [ 1 xf fx ];               % original code
fmin_data = [ ones(size(xf)) xf fx ];    % vectorized code

fv = fx; fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;   tol2 = 2.0*tol1;

d = ones(size(xf)) * -inf;
e = ones(size(xf)) * -inf;
    
% Main loop

options(10) = 0; 
% while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) )         % original version
while ( any(abs(xf-xm) > (tol2 - 0.5*(b-a))) )      % vectorized version

    num = num+1;
    gs = 1;
    I_gs = [1:size(xf,1)];
    
    % Is a parabolic fit possible
    I_poss = find(abs(e) > tol1);
    if ~isempty(I_poss)
        % Yes, so fit parabola
        gs = 0;
        r = (xf-w).*(fx-fv);           % vectorized
        q = (xf-v).*(fx-fw);           % vectorized
        p = (xf-v).*q-(xf-w).*r;       % vectorized
        q = 2.0*(q-r);                 
        
        I_q_pos = find(q > 0.0);
        if ~isempty(I_q_pos)
          p(I_q_pos) = -p(I_q_pos); 
        end
        
        q = abs(q);
        r = e;  
        e = d;

        % Is the parabola acceptable
        I_par = find((abs(p)<abs(0.5.*q.*r)) & (p>q.*(a-xf)) & (p<q.*(b-xf)) );
        I_gs = find((abs(p)>=abs(0.5.*q.*r)) | (p<=q.*(a-xf)) | (p>=q.*(b-xf)) ); 

        if ~isempty(I_par)

            % Yes, parabolic interpolation step
            d(I_par) = p(I_par)./q(I_par);
            x(I_par) = xf(I_par)+d(I_par);
            step = '       parabolic';
     
            % f must not be evaluated too close to ax or bx
            I_close = find(x(I_par)-a(I_par) < tol2(I_par) | ...
                           b(I_par)-x(I_par) < tol2(I_par));
            if ~isempty(I_close)                           
                si = sign(xm(I_par(I_close))-xf(I_par(I_close))) + ...
                         ((xm(I_par(I_close))-xf(I_par(I_close))) == 0);
                d(I_par(I_close)) = tol1(I_par(I_close)).*si;
            end
        end % if ~isempty(I_par)
          
        if ~isempty(I_gs)
            % Not acceptable, must do a golden section step
            gs=1;      
            
        end % if ~isempty(I_gs)
        
    end
    if gs
        % A golden-section step is required
        I_big = find(xf(I_gs) >= xm(I_gs));
        I_small = find(xf(I_gs) < xm(I_gs));

        if ~isempty(I_big)
          e(I_gs(I_big)) = a(I_gs(I_big))-xf(I_gs(I_big));    
        end
        
        if ~isempty(I_small)
          e(I_gs(I_small)) = b(I_gs(I_small))-xf(I_gs(I_small));  
        end
        
        d(I_gs) = c*e(I_gs);
        step = '       golden';
    end

    % The function must not be evaluated too close to xf
    si = sign(d) + (d == 0);
    
    x = xf + si .* max( [abs(d'); tol1'] )';       % vectorized version

    fu = feval(funfcn,x,varargin{:});  
    fmin_data = [num*ones(size(x)) x fu];       % vectorized version

    % Update a, b, v, w, x, xm, tol1, tol2
    I_low = find(fu <= fx);
    I_high = find(fu > fx);

    if ~isempty(I_low)                          % vectorized version 
        I_low2 = find(x(I_low) >= xf(I_low));
        I_not_low2 = find(x(I_low) < xf(I_low));

        if ~isempty(I_low2), 
          a(I_low(I_low2)) = xf(I_low(I_low2));
        end
        if ~isempty(I_not_low2) 
          b(I_low(I_not_low2)) = xf(I_low(I_not_low2)); 
        end
        
        v(I_low) = w(I_low); fv(I_low) = fw(I_low);
        w(I_low) = xf(I_low); fw(I_low) = fx(I_low);
        xf(I_low) = x(I_low); fx(I_low) = fu(I_low); 
    end

    if ~isempty(I_high)                          % vectorized version 
        I_high2 = find(x(I_high) < xf(I_high));
        I_not_high2 = find(x(I_high) >= xf(I_high));

        if ~isempty(I_high2)
          a(I_high(I_high2)) = x(I_high(I_high2));            
        end
        
        if ~isempty(I_not_high2)
        % else,
          b(I_high(I_not_high2)) = x(I_high(I_not_high2)); 
        end
        
        I_1 = find((fu(I_high) <= fw(I_high)) | (w(I_high) == xf(I_high)) );
        if ~isempty(I_1)
            v(I_high(I_1)) = w(I_high(I_1)); 
            fv(I_high(I_1)) = fw(I_high(I_1));
            w(I_high(I_1)) = x(I_high(I_1)); 
            fw(I_high(I_1)) = fu(I_high(I_1));
        end
        
        I_1 = find((fu(I_high) <= fv(I_high)) | (v(I_high) == xf(I_high)) | ...
                   (v(I_high) == w(I_high)) );
        % elseif ( any(fu <= fv) | any(v == xf) | any(v == w) ) 
        if ~isempty(I_1)
            v(I_high(I_1)) = x(I_high(I_1)); 
            fv(I_high(I_1)) = fu(I_high(I_1));
        end
    end

    xm = 0.5*(a+b);
    
    tol1 = seps*abs(xf) + tol/3.0; 
    tol2 = 2.0*tol1;
    
    if num > options(14)
        if options(1) >= 0
        disp('Warning: Maximum number of function evaluations has been');
        disp('         exceeded - increase options(14).')
        options(10)=num;
        return
        end
    end
end
%options(8) = feval(funfcn,x,varargin{:});
options(10)=num;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ELLIPALT function within LOS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fx = ellipalt(los_length, x1, x2, x3, xlos_unit1, xlos_unit2, xlos_unit3);
% Provided a length of the los vector, compute the corresponding tangent altitude
% above WGS 84 ellipsoid.

rt(:,1) = x1 + xlos_unit1 .* los_length;
rt(:,2) = x2 + xlos_unit2 .* los_length;
rt(:,3) = x3 + xlos_unit3 .* los_length;

[lla] = ecef2lla(rt,1);
fx = lla(:,3);

%%%%% END ALGORITHM CODE %%%%%

% end ELLIP_ALT

%%%%% END ALGORITHM CODE %%%%%

% end REMAP

