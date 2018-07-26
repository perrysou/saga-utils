function [dops, t_out, num_sats] = ned2dops(ned,t);

% [dops, t_out, num_sats] = ned2dops(ned,t);
%
% Computes DOP values from North, East, Down LOS vectors.  All LOS data
% that is passed into this function will be used in the DOP computations.
% Masking and data editing should be performed prior to passing the NED
% LOS vectors to NED2DOPS.
%
% Input:
%   ned      - Visible LOS vector in NED coordinates (nx3) unit vectors
%   t        - GPS time corresponding to each row in the LOS vectors, 
%               [GPS_week GPS_sec]  (nx2)
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               valid GPS_sec values are 0-604799
% Output:
%   dops     - dilution of precision for each time (kx5) 
%               where k = number of time steps
%               [GDOP, PDOP, HDOP, VDOP, TDOP]
%   t_out    - output time matrix (kx2) [GPS_week GPS_sec]
%   num_sats - number of satellites visible (kx1)
%
%   Note:  If only 3 satellites are found at a given time, the altitude is
%          assumed fixed and only the HDOP and TDOP are filled with values.
%          If fewer than 3 satellites are found at a given time, the DOPs values
%          are filled with inf for that time.
%
% See also VIS_DATA, LL2DOPS

% Written by: Jimmy LaMance 10/24/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, NORMVECT

% Reference: Parkinson and Spilker (Blue Book) vol. 1, page 413-414. 
% Modified to operate in NED coordinates.

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
dops=[]; t_out=[]; num_sats=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on NED2DOPS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.func_name = 'NED2DOPS';

estruct.variable(1).name = 'ned';
estruct.variable(1).req_dim = [901 3];
estruct.variable(1).var = ned;
  
estruct.variable(2).name = 't';
estruct.variable(2).req_dim = [901 2];
estruct.variable(2).var = t;
estruct.variable(2).type = 'GPS_TIME';
                                
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
% compute the input time in linear format
t_in = t(:,1) * 86400 * 7 + t(:,2);     % convert to seconds past GPS epoch 

% sort the imput in time so all common obs are together
[t1, I] = sort(t_in);
t_in = t_in(I,:);
ned = ned(I,:); 
t = t(I,:);
clear t1 I

% generate the output times 
% find changes in t (different sets of observations)
dt_all = diff(t_in);
dt_all = [1; dt_all];

% search for the time changes
I = find(dt_all ~= 0);
num_times = length(I);

% build the t_out matrix
t_out = t(I,:);

% store some temporary time matrices for later use
t_in_all = t_in;
t_in_temp = t_in(I);

% allocate the space for the dops and num_sats output variable
dops = ones(num_times,5) * inf;     % initialize to infinity
num_sats = zeros(num_times,1);      % initialize to zero

% check that all NED LOS vectors have unit length
if (any(ned' * ned ~= 1))
  ned = normvect(ned);                    % normalize the LOS vectors
end % if

% fill in the time portion of the LOS (soon to be G) matrix
ned(:,4) = ones(size(ned,1),1);  % time portion of the 'Geometry' matrix

% find changes in t (different sets of observations)
dt_all = diff(t_in);
dt_all = [1; dt_all];

% search for the time changes
I = find(dt_all ~= 0);
t_pass = t_in(I);

% put a 1 at all of the time changes
dti_all = dt_all;
dti_all(I) = ones(size(I));

% build the num_sats output (number of satellites for each set/DOP)
svs = [I; size(ned,1)+1];
num_sats_pass = diff(svs);

% put all of the times together
[t_all, I_sort] = sort([t_in_temp; t_pass]);

% put the number of satellites in the same order as the times
num_sats_all = [num_sats; num_sats_pass];
num_sats_all = num_sats_all(I_sort);

% find duplicate times in t_all
dt_all_full = diff(t_all);
dt_all_full = [dt_all_full; length(dt_all_full)];
I_dup = find(dt_all_full == 0);

% build an index to the output num_sats to replace with the vis. satellites
num_sats_pass_index = I_sort(I_dup);

% replace the zero in the num_sats matrix with the visible number of satellites
% from the num_sats_pass matrix.
num_sats(num_sats_pass_index) = num_sats_pass;

% create an index of observation numbers for each set of LOS vectors
% ie (1 1 1 1.... 1 2 2 2 ....2 3 3) etc 
obs_index = cumsum(dti_all);
obs_index_all = obs_index;
ned_all = ned; 

% now build an index the same size as the obs_index that has the number of 
% satellites visible ie (1 4 4 4 4 3 3 3 2 2 5 5 5 5 5...2 2 3 3 3...4 4 4 4) 
% initialize the matrix to be all zeros
num_vis_full = zeros(size(obs_index));

% start but using the dti_all matrix and putting the num_sats_pass on the 
% points where the time changes for a new block of satellites
num_vis_full(I) = dti_all(I) .* num_sats_pass;

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
I_few_sats = find(num_sats < 3);    

% now generate indices into the observations
I_obs_good = find(num_vis_full > 3);
I_obs_3 = find(num_vis_full == 3);
I_obs_few = find(num_vis_full < 3);

% now create observation matrices with 3 and >3 satellites visible
obs_index3 = obs_index(I_obs_3);
obs_index4 = obs_index(I_obs_good);

% fill in the NED vectors for the 3 and >3 satellite cases
ned3 = ned_all(I_obs_3,[1,2,4]);  % remove the height component for 3 sats case
ned4 = ned_all(I_obs_good,:);

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

if any(I_3_sats)
  % generate the index for each of the four columns in ned
  obs_ind1 = (obs_index3 - 1) * 3 + 1;
  obs_ind2 = (obs_index3 - 1) * 3 + 2;
  obs_ind3 = (obs_index3 - 1) * 3 + 3;

  % put together the full index matrix with all four columns
  g_index = [obs_ind1 obs_ind2 obs_ind3];

  % allocate a sparse matrix for G which has G' for each of times
  % one the diagonal
  num_times = length(I_3_sats);
  num_filled_entries = size(ned3,1) * 3;
  total_size_m = size(ned3,1);
  total_size_n = num_times * 3;
  total_size_mn = total_size_m * total_size_n;

  G = spalloc(total_size_m, total_size_n, num_filled_entries);

  % fill in the G matrix
  column_index = reshape(g_index',1,total_size_m*3)';

  row_index1 = [1:total_size_m]' * ones(1,3);
  row_index = reshape(row_index1',total_size_m*3,1);

  % create a linear index so we can do a straight assignment without looping
  g_index = (column_index - 1) * total_size_m + row_index;

  % fill in the diagonal portion of G with the NED LOS vectors
  % (including the time column)
  G(g_index) = ned3';

  % compute the A matrix, the inverse of G' * G
  A = inv(G' * G); 

  % extract the main diagonal from the A-matrix
  sigmas = spdiags(A,0); 

  % the DOPS are a combination of the sigmas, the sigmas are ordered
  % N1, E1, T1, N2, E2, ...T2, N3, ...T3, ...Nm,...Tm
  % so we reshape the sigmas matrix (just to make it easier to work with)
  % into num_times x 3  with the sigmas for each time step filling 
  % in the 3 dimensions
  sig1 = reshape(sigmas,3,num_times)';   

  % now build up the DOPS
  dops(I_3_sats,3) = sqrt(sig1(:,1) + sig1(:,2));       % HDOP
  dops(I_3_sats,5) = sqrt(sig1(:,3));                   % TDOP
  
  % clear out variables so they can be reused for the 4 satellite case
  clear G A g_index obs_ind1 obs_ind2 obs_ind3 num_filled_entries
  clear total_size_m total_size_n total_size_mn 
  clear column_index row_index row_index1
  
end % if any(I_3_sats)  

% if there are not any 4 (or more) satellite cases return
if ~any(I_good)
  return
end % if ~any(I_good)

% generate the index for each of the four columns in ned
obs_ind1 = (obs_index4 - 1) * 4 + 1;
obs_ind2 = (obs_index4 - 1) * 4 + 2;
obs_ind3 = (obs_index4 - 1) * 4 + 3;
obs_ind4 = (obs_index4 - 1) * 4 + 4;

% put together the full index matrix with all four columns
g_index = [obs_ind1 obs_ind2 obs_ind3 obs_ind4];

% allocate a sparse matrix for G which has G' for each of times
% one the diagonal
num_times = length(I_good);
num_filled_entries = size(ned4,1) * 4;
total_size_m = size(ned4,1);
total_size_n = num_times * 4;
total_size_mn = total_size_m * total_size_n;

G = spalloc(total_size_m, total_size_n, num_filled_entries);

% fill in the G matrix
column_index = reshape(g_index',1,total_size_m*4)';

row_index1 = [1:total_size_m]' * ones(1,4);
row_index = reshape(row_index1',total_size_m*4,1);

% create a linear index so we can do a straight assignment without looping
g_index = (column_index - 1) * total_size_m + row_index;

% fill in the diagonal portion of G with the NED LOS vectors
% (including the time column)
G(g_index) = ned4';

% compute the A matrix, the inverse of G' * G
A = inv(G' * G); 

% extract the main diagonal from the A-matrix
sigmas = spdiags(A,0); 

% the DOPS are a combination of the sigmas, the sigmas are ordered
% N1, E1, D1, T1, N2, E2, ...T2, N3, ...T3, ...Nm,...Tm
% so we reshape the sigmas matrix (just to make it easier to work with)
% into num_times x 4  with the sigmas for each time step filling 
% in the 4 dimensions
sig1 = reshape(sigmas,4,num_times)';   

% now build up the DOPS
dops(I_good,1) = sqrt(sig1(:,1) + sig1(:,2) + sig1(:,3) + sig1(:,4));   % GDOP
dops(I_good,2) = sqrt(sig1(:,1) + sig1(:,2) + sig1(:,3));               % PDOP
dops(I_good,3) = sqrt(sig1(:,1) + sig1(:,2));                           % HDOP
dops(I_good,4) = sqrt(sig1(:,3));                                       % VDOP
dops(I_good,5) = sqrt(sig1(:,4));                                       % TDOP

%%%%% END ALGORITHM CODE %%%%%

% end of NED2DOPS
