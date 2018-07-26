function [clk_bias, clk_drift] = clockerr(t_pr,clk_model,seed)

% [clk_bias, clk_drift] = clockerr(t_pr,clk_model,seed);
%
% This function simulates the user clock bias and drift rate of a GPS receiver.
%
% Input:
%   t_pr      - GPS time of pseudoranges (PR) (nx2) [GPS_week GPS_sec]
%                valid GPS_week values are 1-3640 (years 1980-2050)
%                valid GPS_sec values are 0-604799
%   clk_model - input for the user clock model (optional)
%                1x2 matrix in the form [S_b S_f] (1x2)
%                where the S_b and S_f are white noise spectral 
%                densities in s and s/s. These values can easily 
%                be related to Allan variances (see manual for details).
%                Default is the based on a crystal oscillator with
%                values of [4e-19 1.58e-18]
%   seed      - seed value for random number generator (optional)
%                Default value is 0.
% Output:
%   clk_bias  - User clock bias (m) (nx1)
%   clk_drift - User clock drift rate (m/s) (nx1)
%
% See also PSEUDO_R, SA_EPS, SA_CLOCK, TROPDLAY

% Written by: Jimmy LaMance 1/15/97
% Copyright (c) 1998 by Constell, Inc.

% Reference: Global Positioning System: Theory and Applications
% Volume 1, Parkinson and Spilker, pages 417-418.

% functions called: ERR_CHK

% WGS-84 constants
RE = 6378137;                  % WGS-84 value in meters
LIGHT_SPEED = 299792458;       % WGS-84 value in m / s

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
clk_bias=[]; clk_drift=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on CLOCKERR for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 2
  clk_model = [4e-19 1.58e-18];
end

if nargin < 3
  seed = 0;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'CLOCKERR';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_pr;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'clk_model';
  estruct.variable(2).req_dim = [1 2];
  estruct.variable(2).var = clk_model;
  
  estruct.variable(3).name = 'seed';
  estruct.variable(3).req_dim = [1 1];
  estruct.variable(3).var = seed;
  
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
                                                                              
% allocate the output parameters.  this has valid and invalid times.
clk_bias = ones(size(t_pr,1),1) * inf;
clk_drift = ones(size(t_pr,1),1) * inf;

% seed the random number generator
randn('seed',seed);

% load clock model data into local variables
Sb = clk_model(1);
Sf = clk_model(2);

% get everything in seconds past a common minimum week
t = gpst2sec(t_pr);   

% Set a sort flag and check to see if the data need to be time sorted
sort_flag = 0;
if any(diff(t) < 0)
  [t, I_sort] = sort(t);
  sort_flag = 1;
end % if any(diff(t) < 0)

% sort the data so that time is always increasing
% keep track of the index used to sort the times, so we can return everything
% in the same order it came in
[tu, I_t1] = unique(t);

% start the time matrix t at 0
tu = tu - tu(1);

% Sort out the unique times and create indices to all of the input times.
t1_all_index = zeros(size(t));

num_t1 = length(I_t1);

% Remove the last intersection point from both index data sets
I_t1 = I_t1(1:end-1);

% Insert a 1 at all of the intersection time changes.
t1_all_index(I_t1+1) = ones(size(I_t1));

% Put a 1 at the starting index
t1_all_index(1) = 1;

% Use the cumsum function to create the indices to the original 
% times.  This will results in matrices like [1 1 1 1 2 2 2 2... n n n].
t1_all_index = cumsum(t1_all_index);

% make sure that there is more than one time input
if length(t) == 1    % if there is only 1 time to process
  % return a warning message and no output 
  clk_bias = 0;
  clk_drift = 0;
  fprintf('There was only a single time sent to CLOCKERR.\n');
  fprintf('The function is designed for multiple time inputs to generate\n');
  fprintf('correlated clock errors.  See help on CLOCKERR for details.\n');  
  fprintf('Returning to calling function with blank output.\n');
  return
end % if length(t) <= 1 
  
% compute the delta-t for use in the algorithm
dt = diff(tu);

% verify that the dt value is never zero
if any(dt == 0)
  % return a warning message and no output 
  clk_bias = 0;
  clk_drift = 0;
  fprintf('There was a problem creating unique times from, t_pr,');
  fprintf('for CLOCKERR.\n');
  fprintf('The function is designed to use monotomically increasing time inputs to\n');
  fprintf('generate correlated clock errors.  See help on CLOCKERR for details.\n');  
  fprintf('Contact technical support.\n');
  fprintf('Returning to calling function with empty output.\n');
  return
end % if any(dt==0) 

% compute the w matrix, w will be (n x 2)
num_time_steps = length(tu);

sig_bias = sqrt(Sb*dt(1));      % 1-sigma value of clock bias noise (sec)
sig_drift = sqrt(Sf*dt(1));     % 1-sigma value of clock drift noise (sec/sec)

% Build up a transformed white noise model (as per Parkinson, et al)

wn_bias = randn(num_time_steps,1)*sig_bias;     % vector of bias noise
wn_drift = randn(num_time_steps,1)*sig_drift;   % vector of drift noise

wnc = [wn_bias, wn_drift]*[1 dt(1)/2; 0 1];     % transformed noise matrix

% If the time step is constant, and ltitr exists, use the vectorized version
% of the code

% Find out if the time intervals are uniform

diff_dt = diff(dt);     % difference in time step size
if any(abs(diff_dt - diff_dt(1)) > 1e-4*diff_dt(1)),
    const_step = 0;     % non-constant time step size
else
    const_step = 1;     % constant time step size
end;

% If step size is constant and the function ltitr exists on the current 
% machine, vectorize the clock noise process

if ((const_step == 1) & (exist('ltitr') == 5)), % vectorize code

    x0 = [randn(1,1)*5e-5,randn(1,1)*5e-8]';

    A = [1 dt(1); 0 1];
    B = [1 0; 0 1];

    y = ltitr(A,B,wnc,x0);

    clk_bias = y(:,1)*LIGHT_SPEED;
    clk_drift = y(:,2)*LIGHT_SPEED;

else,                                           % non-vectorized code

    % seed the process with the first value
    clk_bias(1) = randn(1,1) * 5e-5;  % 50 msec 1-sigma of initial error
    clk_drift(1) = randn(1,1) * 5e-8; % 50 ns/s 1-sigma of initial error 

    % loop over the time steps
    for i = 2:num_time_steps,

      count = i - 1;
      
      % update clock drift and bias for this time
      clk_bias(i) = clk_bias(count) + dt(count) .* ...
                            clk_drift(count) + wnc(count,1);
      clk_drift(i) = clk_drift(count) + wnc(count,2);
    end % for i = 2:num_time_steps 
           
    % convert to m and m / s, and get back in the same order as inputs
    clk_bias = clk_bias * LIGHT_SPEED;
    clk_drift = clk_drift * LIGHT_SPEED;

end;

% Expand the clock bias and drift to be the same size and ordering as the input
% times.
clk_bias = clk_bias(t1_all_index);
clk_drift = clk_drift(t1_all_index);

% Unsort the data if the sort_flag has been set
if sort_flag == 1
  clk_bias = clk_bias(I_sort);
  clk_drift = clk_drift(I_sort);
end % if sort_flag == 1       

%%%%% END ALGORITHM CODE %%%%%

% end of CLOCKERR
