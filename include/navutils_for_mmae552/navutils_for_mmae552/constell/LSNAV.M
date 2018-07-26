function [t_nav,x_nav,num_sats,nav_index,v_nav] = ...
                    lsnav(t_pr,pr,gps_orb,init_pos,doppler,gps_vel,init_vel);

% [t_nav,x_nav,num_sats,nav_index] = lsnav(t_pr,pr,gps_orb,init_pos); 
%      or for a position and velocity solution 
% [t_nav,x_nav,num_sats,nav_index,v_nav] = ...
%         lsnav(t_pr,pr,gps_orb,init_pos,doppler,gps_vel,init_vel);
%
% Computes a least squares PR navigation (position) solution or a
% least squares position and velocity solution if the Doppler data is provided.
%
% Input:
%   t_pr      - time associated with the PR and orbits [GPS_week GPS_sec] (kx2)
%                valid GPS_week values are 1-3640 (years 1980-2050)
%                valid GPS_sec values are 0-604799
%   pr        - pseudo-range (PR) corresponding to t_pr time (kx1) (m) 
%   gps_orb   - GPS/GLONASS satellite orbits corresponding to pr 
%                (kx4) (m) ECEF [sv_num x y z]
%   init_pos  - initial estimate of the user position and clock bias
%                (1x4) (m) ECEF [x y z clk_bias]
%   doppler   - Doppler measurements corresponding to t_pr time (kx1) (m)
%                (optional).  No velocity solution is returned if not provided.
%   gps_vel   - GPS/GLONASS satellite velocities corresponding to pr and 
%                Doppler measurements (kx3) (m/s) ECEF [vx vy vz] 
%                (required when using Doppler data)
%   init_vel  - initial estimate of the user velocity and clock drift
%                (1x4) (m/s) ECEF [vx vy vz clk_drift]
%                (required when using Doppler data)
%
% Output:
%   t_nav     - time of navigation solutions (nx2) [GPS_week GPS_sec]
%   x_nav     - navigation position solution for each time (nx4) (m) ECEF
%                [x y z clk_bias]
%   num_sats  - number of satellites used in the nav computation (nx1)
%                n is the number of resulting navigation solutions
%   nav_index - index that relates the t_nav matrix to the t_pr matrix (kx1)
%                ie [1 1 1...2 2 2...3 3 3 3 3...] where all of the 1s refer
%                to the first t_nav/t_pr, the 2s the the next, etc.
%   v_nav     - navigation velocity solution for each time (nx4) (m/s) ECEF
%                [vx vy vz clk_drift]
%
% Note: The PR modeling contained in LSNAV does not model the satellite motion,
%       satellite clock, line bias, or relativity.  These effects should be 
%       removed from the PR measurements before calling LSNAV.
% Note: Times without 3 or more visibile satellites will be filled with the 
%       initial position and velocity estimates provided in init_pos and
%       init_vel.
% Note: For initial altitude estimates of greater than 10 km, the vehicle 
%       trajectory is assumed to be a satellite and a fixed altitude solution
%       with only 3 satellites is not computed.  All other trajectories will
%       compute an alititude hold navigation solution when only 3 satellites
%       are available.
%
% See also PROPGEPH, PSEUDO_R, DOPPLER

% Written by: Jimmy LaMance 2/28/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPST2SEC, ECEF2LLA, ECEF2NED, NED2ECEF,
%                   NORMVECT

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
t_nav=[]; x_nav=[]; num_sats=[]; nav_index=[]; v_nav=[];

% Check the number of input arguments and issues a message if invalid
if nargin ~= 4 & nargin ~= 7
  fprintf('Invalid number of input variable to LSNAV.  4 or 7 inputs required.\n');
  fprintf('Type "help lsnav" for details.\n\n');
  return
end % if nargin ~= 4 & nargin ~= 7

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'LSNAV';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_pr;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'pr';
  estruct.variable(2).req_dim = [901 1];
  estruct.variable(2).var = pr;
  
  estruct.variable(3).name = 'gps_orb';
  estruct.variable(3).req_dim = [901 4];
  estruct.variable(3).var = gps_orb;

  estruct.variable(4).name = 'init_pos';
  estruct.variable(4).req_dim = [1 4];
  estruct.variable(4).var = init_pos;
  
  if nargin == 7
    estruct.variable(5).name = 'doppler';
    estruct.variable(5).req_dim = [901 1];
    estruct.variable(5).var = doppler;
  
    estruct.variable(6).name = 'gps_vel';
    estruct.variable(6).req_dim = [901 3];
    estruct.variable(6).var = gps_vel;
  
    estruct.variable(7).name = 'init_vel';
    estruct.variable(7).req_dim = [1 4];
    estruct.variable(7).var = init_vel;
  end % if nargin == 7
  
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

% compute the input time in linear format, seconds past GPS epoch
[t_in] = gpst2sec(t_pr);     

% sort the imput in time so all common obs are together
[t_in, I] = sort(t_in);
t_pr = t_pr(I,:);

% resort the remained of the inputs also
pr = pr(I);
gps_orb = gps_orb(I,:);

if nargin == 7
  doppler = doppler(I);
  gps_vel = gps_vel(I,:);
end % if nargin == 7
  
clear t1 I

% generate the output times 
% find changes in t (different sets of observations)
dt_all = diff(t_in);
dt_all = [1; dt_all];

% search for the time changes
I = find(dt_all ~= 0);
num_times = length(I);

% build the t_out matrix
t_nav = t_pr(I,:);

% store some temporary time matrices for later use
t_in_all = t_in;
t_in_temp = t_in(I);

% allocate the space for the dops and num_sats output variable
x_nav = ones(num_times,1) * init_pos;

if nargin == 7
  v_nav = ones(num_times,1) * init_vel;
else
  v_nav = []; 
end % if nargin == 7

% put a 1 at all of the time changes
dti_all = dt_all;
dti_all(I) = ones(size(I));

% build the num_sats output (number of satellites for each set/DOP)
svs = [I; size(t_pr,1)+1];
num_sats = diff(svs);

% create an index of observation numbers for each set of LOS vectors
% ie (1 1 1 1.... 1 2 2 2 ....2 3 3) etc 
obs_index = cumsum(dti_all);
obs_index_all = obs_index;

% now build an index the same size as the obs_index that has the number of 
% satellites visible ie (4 4 4 4 3 3 3 2 2 5 5 5 5 5...2 2 3 3 3...4 4 4 4) 
% initialize the matrix to be all zeros
num_vis_full = zeros(size(obs_index));

% start but using the dti_all matrix and putting the num_sats_pass on the 
% points where the time changes for a new block of satellites
num_vis_full(I) = dti_all(I) .* num_sats;

% now modify the transition points such that the n+1 transition value
% is equal to the current n+1 transition less the n transition
length_I = length(I);
num_vis_full(I(2:length_I)) = ...
      num_vis_full(I(2:length_I)) - num_vis_full(I(1:length_I-1));     
      
% finish off the process using the cumcum function to fill in all of the zeros
% to the appropriate value with the number of satellites visible
num_vis_full = cumsum(num_vis_full);      

% find time steps when there are only 3 satellites or less than 3 satellites 
% for >= 4 satellites we will do a inverse of the A matrix using sparse 
% matrix functions
% for == 3 satellites we will loop over each A matrix for the case and do
% a full matrix pinv
% for < 3, just return the defualt values of inf for the dops, but leave
% the num_sats variable filled correctly

% start with the index into the DOP solutions
I_good = find(num_sats > 3);
I_3_sats = find(num_sats == 3);    

% now generate indices into the observations
I_obs_good = find(num_vis_full > 3);       % 4 or more satellites
I_obs_3 = find(num_vis_full == 3);         % 3 satellites
I_obs_nav = find(num_vis_full >= 3);       % all satellites used in the nav

% now create observation matrices with 3 and >3 satellites visible
obs_index3 = obs_index(I_obs_3);
obs_index4 = obs_index(I_obs_good);

nav_index = obs_index(I_obs_nav);

% save the obs_index values into an obs_index_nav3 and obs_index_nav4 matrices
% because the obs_index? matrices will be modified below
obs_index_nav_3 = obs_index3;
obs_index_nav_4 = obs_index4;

% reduce the data to 3 and >3 version of the data for storage for later
% use.  This enables the most efficient use of the sparse matrix allocation
pr3 = pr(I_obs_3);
pr4 = pr(I_obs_good);

gps_orb3 = gps_orb(I_obs_3,:);
gps_orb4 = gps_orb(I_obs_good,:);

if nargin == 7
  doppler3 = doppler(I_obs_3);
  doppler4 = doppler(I_obs_good);

  gps_vel3 = gps_vel(I_obs_3,:);
  gps_vel4 = gps_vel(I_obs_good,:);

end % if nargin == 7
  
% renumber the obs_index for the two cases such that the index is sequential
% starting at 1 (no skips)
I_diff_3 = [1; (diff(obs_index3))];
I_diff_4 = [1; (diff(obs_index4))];

% sort out any diffs that are greater than 1 and set them to 1
I_big_diff3 = find(I_diff_3 > 1);
I_big_diff4 = find(I_diff_4 > 1);

if ~isempty(I_big_diff3)
  I_diff_3(I_big_diff3) = ones(size(I_diff_3(I_big_diff3)));
end % if ~isempty(I_big_diff3) 

if ~isempty(I_big_diff4)
  I_diff_4(I_big_diff4) = ones(size(I_diff_4(I_big_diff4)));
end % if ~isempty(I_big_diff4) 

% now build up the new obs_index
obs_index3 = cumsum(I_diff_3);
obs_index4 = cumsum(I_diff_4);

% initialize variables for the least squares
max_iter = 10;           % maximum number of LS iterations

converge_flag = 0;       % flag for overall convergence
x_converge_flag = 0;     % flag for position convergence
v_converge_flag = 0;     % flag for velocity convergence

dx_converge = .01;       % position convergence criterion
dv_converge = .01;        % velocity converge criterion
i_count = 1;             % initialize the nav solution counter
iter = 0;      % iteration count
dx = 100000;   % initialization of convergence

% check to see if this is a satellite soltuion (height > 10 km)
test_lla = ecef2lla(x_nav(1,1:3));

% if it is a satellite, don't try to do any 2-D navigation solutions
if test_lla(3) > 6378137 + 1e5
  orbital_flag = 1;
else
  orbital_flag = 0;
end % if

if any(I_3_sats) & orbital_flag ~= 1
  % generate the index for each of the four columns in ned
  obs_ind1 = (obs_index3 - 1) * 3 + 1;
  obs_ind2 = (obs_index3 - 1) * 3 + 2;
  obs_ind3 = (obs_index3 - 1) * 3 + 3;

  % put together the full index matrix with all four columns
  g_index = [obs_ind1 obs_ind2 obs_ind3];

  % allocate a sparse matrix for G which has G' for each of times
  % one the diagonal
  num_obs = length(obs_index3);
  num_times = length(I_3_sats);
  num_filled_entries = num_obs * 3;
  total_size_m = num_obs;
  total_size_n = num_times * 3;
  total_size_mn = total_size_m * total_size_n;

  G = spalloc(total_size_m, total_size_n, num_filled_entries);

  % fill in the G matrix
  column_index = reshape(g_index',1,total_size_m*3)';

  row_index1 = [1:total_size_m]' * ones(1,3);
  row_index = reshape(row_index1',total_size_m*3,1);

  % create a linear index so we can do a straight assignment without looping
  g_index = (column_index - 1) * total_size_m + row_index;

  % set the iteration counter and the converge flag
  iter = 0;
  x_converge_flag = 0;                               
  v_converge_flag = 0;                               
    
  % loop until the maximum number of iteration has been exceeded or all
  % of the navigation solutions have converged
  while iter < max_iter & converge_flag == 0

    % compute the LOS vector to the GPS satellites for the trajectory
    los_ecef = x_nav(obs_index_nav_3,1:3) - gps_orb3(:,2:4);    
    
    % normalize the LOS vector
    [los_ecef, r_mag] = normvect(los_ecef);

    % convert ECEF LOS vectors to NED for fixed altitude computations
    [los_ned] = ecef2ned(los_ecef, ecef2lla(x_nav(obs_index_nav_3,1:3)));

    % fill in the time portion (and zero out the height component)
    los_ned(:,3) = ones(size(los_ned,1),1);  % time portion of the 'Geometry' matrix
      
    % fill in the diagonal portion of G with the NED LOS vectors
    % (including the time column)
    G(g_index) = los_ned';
    
    % compute the computed PR measurement with the clock term
    comp_pr = r_mag + x_nav(obs_index_nav_3,4);  % clock term

    % compute residuals
    z = pr3 - comp_pr;

    % compute the temporary value inv(G'G)G'.  this will be used again 
    % in the Doppler processing
    GtG_temp = inv(G' * G) * G';
    
    % do least squares inv(G' * G) * G' * z
    delta_x_4(:,[1,2,4]) = reshape((GtG_temp * z)',3,num_times)';   
      
    % add the height component back into delta_x_nav
    delta_x_4(:,3) = zeros(size(delta_x_4,1),1);

    % rotate the delta_x_nav vector back to ECEF and add the clock term
    delta_x_nav = ned2ecef(delta_x_4(:,1:3), ecef2lla(x_nav(I_3_sats,1:3)));
    delta_x_nav(:,4) = delta_x_4(:,4);   % clock term
      
    % convergence calculation (including the clock)
    [dx_vect, dx_mag] = normvect(delta_x_nav);    % delta_x_nav is (1x4)
    
    % check for convergence
    I = find(dx_mag > dx_converge);
    if ~any(I)       
      x_converge_flag = 1;
    end % if ~any(I)
  
    clear I
    % update the navigation solution (position & clock)
    x_nav(I_3_sats,:) = x_nav(I_3_sats,:) + delta_x_nav;

    % if there is doppler data, also do a velocity solution
    if nargin == 7
      
      % compute the relative velocity vector
      vel_ecef = v_nav(obs_index_nav_3,1:3) - gps_vel3;
      
      % compute the computed Doppler measurement with the clock term
      comp_dop = dot(los_ecef',vel_ecef')' + v_nav(obs_index_nav_3,4);  
     
      % compute the measurement residual 
      z = doppler3 - comp_dop;

      delta_v_4(:,[1,2,4]) = reshape((GtG_temp * z)',3,num_times)';   
      
      % add the height component back into delta_x_nav
      delta_v_4(:,3) = zeros(size(delta_v_4,1),1);
      
      % rotate the delta_x_nav vector back to ECEF and add the clock term
      delta_v_nav = ned2ecef(delta_v_4(:,1:3), ecef2lla(x_nav(I_3_sats,1:3)));
      delta_v_nav(:,4) = delta_v_4(:,4);   % clock term
      
      % convergence calculation (including the clock)
      [dv_vect, dv_mag] = normvect(delta_v_nav);    % delta_v_nav is (1x4)
     
      % check for convergence
      I = find(dv_mag > dv_converge);
      if ~any(I)       
        v_converge_flag = 1;
      end % if ~any(I)
  
      clear I
      % update the navigation solution (position & clock)
      v_nav(I_3_sats,:) = v_nav(I_3_sats,:) + delta_v_nav;

    end % if nargin == 7
    
    % set the overall loop converge flag
    if nargin == 7
      if x_converge_flag == 1 & v_converge_flag == 1       
        converge_flag = 1;
      end % if x_converge_flag == 1 & v_converge_flag == 1
    
    else
      if x_converge_flag == 1
         converge_flag = 1;
      end % if x_converge_flag == 1
    end % if nargin == 7
         
    iter = iter + 1;    % increase iteration count
    
    if iter == max_iter
      fprintf('Warning from LSNAV.  Maximum number of iterations exceeded\n');
      fprintf('3 satellite cases.  Results may be in error.\n');
        
      % check to see if convergence failed or FAILED and fix x_nav accordingly
      if dx_mag > 100 * dx_converge
        x_nav(I_3_sats,:) = ones(length(I_3_sats),4) * inf;
        fprintf('Warning from LSNAV.  3 satellite solutions not converged to\n');
        fprintf('within 100 times the convergence the convergence criterion.\n');
        fprintf('All 3 satellite cases will be filled with inf.\n\n');
      end % if dx_mag > 100 * dx_converge3
        
    end % if
    
    end % while iter < max_iter & converge_flag == 0 
    
    % clear out this A and G, will use the same variables later (sparse)
    clear G delta_x_nav delta_x_nav4 delta_x_nav_ned delta_x_nav_ned4 I_i

end % if any(I_3_sats) & orbital_flag ~= 1 

% if there are not any 4 (or more) satellite cases return
if ~any(I_good)
  return
end % if ~any(I_good)

% allocate a sparse matrix for storing the G matrix
% generate the index for each of the four columns in G
obs_ind1 = (obs_index4 - 1) * 4 + 1;
obs_ind2 = (obs_index4 - 1) * 4 + 2;
obs_ind3 = (obs_index4 - 1) * 4 + 3;
obs_ind4 = (obs_index4 - 1) * 4 + 4;

% put together the full index matrix with all four columns
g_index = [obs_ind1 obs_ind2 obs_ind3 obs_ind4];

% allocate a sparse matrix for G
% one the diagonal 
num_obs = length(obs_index4);
num_times = length(I_good);
num_filled_entries = num_obs * 4;
total_size_m = num_obs;
total_size_n = num_times * 4;
total_size_mn = total_size_m * total_size_n;

% allocate the sparse matrix
G = spalloc(total_size_m, total_size_n, num_filled_entries);

% fill in the G matrix
column_index = reshape(g_index',1,total_size_m*4)';

row_index1 = [1:total_size_m]' * ones(1,4);
row_index = reshape(row_index1',total_size_m*4,1);

% create a linear index so we can do a straight assignment without looping
g_index = (column_index - 1) * total_size_m + row_index;
           
% loop while any of the solutions have not converged and have not exceeded
% the maximum number of iterations                     
iter = 0;
converge_flag = 0;       % flag for overall convergence
x_converge_flag = 0;     % flag for position convergence
v_converge_flag = 0;     % flag for velocity convergence

while iter < max_iter & converge_flag == 0
  % compute line of sight vectors
  los_ecef = x_nav(obs_index_nav_4,1:3) - gps_orb4(:,2:4);    

  % normalize the LOS vector
  [los_norm, rmag] = normvect(los_ecef);                                    
  
  % fill in the diagonal portion of G with the NED LOS vectors
  % (including the time column)
  % fill in the time portion of the LOS (soon to be G) matrix
  los_norm(:,4) = ones(size(los_ecef,1),1);  

  G(g_index) = los_norm';

  % compute the computed PR measurement with the clock term
  comp_pr = rmag + x_nav(obs_index_nav_4,4);  % clock term
    
  % compute residuals
  z = pr4 - comp_pr;
                   
  % compute the temporary value inv(G'G)G'.  this will be used again 
  % in the Doppler processing
  GtG_temp = inv(G' * G) * G';

  % least squares navigation solution
  delta_x_nav = reshape((GtG_temp * z),4,num_times)';   

  % convergence calculation (including the clock)
  [dx_vect, dx_mag] = normvect(delta_x_nav);    % delta_x_nav is (1x4)
    
  % check for convergence
  I = find(dx_mag > dx_converge);
  if ~any(I)       
    converge_flag = 1;
  end % if ~any(I)

  clear I

  % update the navigation solution (position & clock)
  x_nav(I_good,:) = x_nav(I_good,:) + delta_x_nav;

  % if there is doppler data, also do a velocity solution
  if nargin == 7
      
    % compute the relative velocity vector
    vel_ecef = v_nav(obs_index_nav_4,1:3) - gps_vel4;
      
    % compute the computed PR measurement with the clock term
    comp_dop = dot(los_norm(:,1:3)',vel_ecef')' + v_nav(obs_index_nav_4,4);  
     
    % compute the measurement residual
    z = doppler4 - comp_dop; 

    delta_v_nav = reshape((GtG_temp * z)',4,num_times)';   
      
    % convergence calculation (including the clock)
    [dv_vect, dv_mag] = normvect(delta_v_nav);    % delta_v_nav is (1x4)
     
    % check for convergence
    I = find(dv_mag > dv_converge);
    if ~any(I)       
      v_converge_flag = 1;
    end % if ~any(I)
  
    clear I
    % update the navigation solution (position & clock)
    v_nav(I_good,:) = v_nav(I_good,:) + delta_v_nav;
      
  end % if nargin == 7
    
  % set the overall loop converge flag
  if nargin == 7
    if x_converge_flag == 1 & v_converge_flag == 1       
      converge_flag = 1;
    end % if x_converge_flag == 1 & v_converge_flag == 1
    
  else
    if x_converge_flag == 1
       converge_flag = 1;
    end % if x_converge_flag == 1
  end % if nargin == 7

  iter = iter + 1;    % increase iteration count
    
  if iter == max_iter
    fprintf('Warning from LSNAV.  Maximum number of iterations exceeded\n');
    fprintf('in 4 satellite case.  Results may be in error.\n');
  end % if
    
end % while iter < max_iter & converge_flag == 0 

%%%%% END ALGORITHM CODE %%%%%

% end of LSNAV
