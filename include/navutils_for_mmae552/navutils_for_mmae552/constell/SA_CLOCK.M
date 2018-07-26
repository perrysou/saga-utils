function [PR_error, PRR_error] = sa_clock(t_pr,prn,sa_model_data,seed)

% [PR_error, PRR_error] = sa_clock(t_pr,prn,sa_model_data,seed);
%
% Simulate the SA clock dither contribution to the pseudo-range (PR), 
% accumulated acrrier phase (CPH) and Doppler measurement models using 
% a 2-nd order Gauss-Markov process.
%
% Input:
%   t_pr           - GPS time of pseudoranges (PR) or carrier phase (CPH)
%                     measurements, (nx2) [GPS_week GPS_sec]
%                     valid GPS_week values are 1-3640 (years 1980-2050)
%                     valid GPS_sec values are 0-604799
%   prn            - GPS satellite numbers corresponding to t_pr (nx1) 
%   sa_model_data  - input for the dither model (optional)
%                     1x3 (or 3x1) matrix in the form [sigma_pr sigma_prr tau]
%                     where the sigmas on the PR and PR-rate are in
%                     m and m/s and tau is the decorrelation time
%                     in seconds.  Default is the RTCA proposed
%                     values of [23 .28 118]
%   seed           - seed value for random number generator (optional)
%                     Default value is 0.
% Output:
%   PR_error       - Pseudo-range error (nx1) (m)
%                     n = num_time_steps x num_sats 
%   PRR_error      - Pseudo-range rate error (nx1) (m/s)
%
% Note: For the 2-nd order Gauss-Markov process the input PRs should be
% evenly spaced in time.  If they are not a warning message will be issued
% and the inputs will be assumed to be evenly spaced.  This will result in
% PR errors that are skewed based on the time tags.
%
% See also PSEUDO_R, SA_EPS, CLOCKERR, TROPDLAY

% Written by: Jimmy LaMance 1/14/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPST2SEC

% Reference: Global Positioning System: Theory and Applications
% Volume 1, Parkinson and Spilker, pages 605-608.

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
PR_error=[]; PRR_error=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,4,nargin);
if ~isempty(msg)
  fprintf('%s  See help on SA_CLOCK for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check that the size of sa_model_data is 1x3, if provided
if nargin >= 3            
  use_RTCA = 0; % use the user provided inputs, (not RTCA standard)
else               % if the sa_model_data is not provided
  use_RTCA = 1;    % use the RTCA standard
  sa_model_data = [23 .28 118];
end % if nargin >= 3

% check that the size of seed is 1x1, if provided
if nargin < 4            
  seed = 0;      % use the default value of 0
end % if nargin < 4

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'SA_CLOCK';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_pr;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'prn';
  estruct.variable(2).req_dim = [901 1];
  estruct.variable(2).var = prn;
  
  estruct.variable(3).name = 'sa_model_data';
  estruct.variable(3).req_dim = [1 3];
  estruct.variable(3).var = sa_model_data;

  estruct.variable(4).name = 'seed';
  estruct.variable(4).req_dim = [1 1];
  estruct.variable(4).var = seed;

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

% sort the time in increasing order, start by getting time into a linear format
tl = gpst2sec(t_pr);
[tl, I_sort] = sort(tl);

% reorder the input time and prn matrices to be in time order
t_pr = t_pr(I_sort,:);
prn = prn(I_sort);

% allocate the output data
PR_error = zeros(size(t_pr,1),1);
PRR_error = zeros(size(t_pr,1),1);

% set up constants fo use in computing the errors
if (use_RTCA == 1)
  % RTCA proposed constants
  beta = 1 / sqrt(2);
  omega_0 = .012;        % rad / sec
  c_squared = .002528;   % m^2    
  sa_model_data = [23 .28 118];
else
  % compute constants based on inputs of sigma_pr, sigma_cph
  % and the time constant tau from the sa_model_data variable
  beta = sa_model_data(1) / (sa_model_data(2) * sa_model_data(3));
  omega_0 = 1 / (beta * sa_model_data(3));
  c_squared = sa_model_data(2)^2  * 4 * beta * omega_0;
end % if

% compute the damped frequency
omega_1 = omega_0 * sqrt(1 - beta^2);

% determine the total number of satellites
% start by sorting the prn numbers and looking for changes
prn_sort = sort(prn);
prn_change = find(diff(prn_sort) ~= 0);

% create a matrix that has sorted and reduced prns [1 2 4 5 6 8 ... 28]
prn_list = [prn_sort(prn_change); prn_sort(length(prn_sort))];

% compute the total number of prns 
num_prns = length(prn_list);

% seed the random number generator
randn('seed',seed);

% Note on the upcoming for loop.
% Compute the S/A noise for each satellite.
% Each satellite has correlated noise, but is decorrelated
% between satellites, so we have to do this 1 satellite at a time.

% loop over the satellite numbers
for i = 1:num_prns
  % find all of the data for this satellite
  I = find(prn == prn_list(i));
  if ~any(I)
    fprintf('Warning in processing PR measurements in SA_CLOCK.\n')
    fprintf('A satellite has been lost in the processing.\n');
    fprintf('Output may be in error.\n\n');
  end % if any(I)

  % make sure that there is more than 1 observation before getting into the next
  % set of algorithms
            
  if length(I) == 1
    % pick a random number (with the right sigma) for the SA contribution
    PR_error(I) = randn(1,1) * sa_model_data(1);
    PRR_error(I) = randn(1,1) * sa_model_data(2);
  
  else % if there are more than 1 observation for this satellite

    % compute the time step between PR measurements.  This must be a constant
    % time step for the algorithm to work correctly.  Use the minimum dt for
    % this satellite.  This is done using the linear time matrix computed above.
    dt = min(diff(tl(I)));   % min time difference for this sat. in GPS sec

    % compute the state transition matrix
    phi(1,1) = exp(-beta * omega_0 * dt) * ...
               (cos(omega_1 * dt) + ...
               beta * (omega_0 / omega_1) * sin(omega_1 * dt));
    phi(1,2) = (1 / omega_1) * exp(-beta * omega_0 * dt) * ...
                sin(omega_1 * dt);
    phi(2,1) = -omega_0^2 * phi(1,2);
    phi(2,2) = exp(-beta * omega_0 * dt) * ...
               (cos(omega_1 * dt) - ...
               beta * (omega_0 / omega_1) * sin(omega_1 * dt));

    % Compute covariance matrix. First create variables for combinations of
    % parameters that are used more than once in the computations
    const_1 = c_squared / (4 * beta * omega_0^3);
    const_2 = exp(-2 * beta * omega_0 * dt);
    ratio2 = omega_0^2 / omega_1^2;
    rate = 2 * omega_1 * dt;

    % build up the Q matrix
    Q(1,1) = const_1 * (1 - ratio2 * const_2 * ...
             (1 - beta^2 * cos(rate) + ...
             beta * (omega_1 / omega_0) * sin(rate)));
    Q(1,2) = (c_squared / (4 * omega_1^2)) * ...
             (const_2 * (1 - cos(rate)));
    Q(2,2) = const_1 * omega_0 ^ 2 * (1 - ratio2 * const_2 * ...
             (1 - beta^2 * cos(rate) - ...
             beta * (omega_1 / omega_0) * sin(rate)));

    % form UD decomposition of Q
      u = zeros(2);
    u(1,1) = sqrt(Q(1,1) - Q(1,2)^2 / Q(2,2));
    u(1,2) = Q(1,2) / sqrt(Q(2,2));
    u(2,2) = sqrt(Q(2,2));

    % compute the number of PR measurements we're working with for this prn
    num_prs = length(I);
   
    % generate zero mean, unit variance, Gaussian noise matrices
    % for the PR and CPH measurements 
    rand_noise = randn(2,num_prs);
    PR_noise = rand_noise(1,:)';
    PRR_noise = rand_noise(2,:)';

    % pick a random number (with the right sigma) to start the process.
    % This will start the process at a non-zero quantity
    PR_error(I(1)) = randn(1,1) * sa_model_data(1);
    PRR_error(I(1)) = randn(1,1) * sa_model_data(2);

        if (exist('ltitr') == 5),   % vectorize code if ltitr function is available

            x0 = [PR_error(I(1)),PRR_error(I(1))]';  % initial noise value

            A = phi;
            B = u;
            Noise = [PR_noise,PRR_noise];           % array of noise inputs to sa model

            y = ltitr(A,B,Noise,x0);                % generate the sa_clock data

            PR_error(I,:) = y(:,1);
            PRR_error(I,:) = y(:,2);

        else,                               % non-vectorized version

            % loop over these PRs

            for i = 2:num_prs
                % compute indices to the PR_error and markov process matrices 
                PR_error(I(i),:) = phi(1,1) * PR_error(I(i-1),:) + ...
                                   phi(1,2) * PRR_error(I(i-1),:) + ...
                                   u(1,1) * PR_noise(i,:) + ...
                                   u(1,2) * PRR_noise(i,:);
                PRR_error(I(i),:) = phi(2,1) * PR_error(I(i-1),:) + ...
                                    phi(2,2) * PRR_error(I(i-1),:) + ...
                                    u(2,2) * PRR_noise(i,:);
            end % for

        end;   % end of if (exist('ltitr'))

  end  % if length(I) == 1
end % for i = 1:length(prn)

%%%%% END ALGORITHM CODE %%%%%

% end SA_CLOCK

