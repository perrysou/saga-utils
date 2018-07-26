function [pr_sa_eps, sa_eps_err] = sa_eps(pr, sigma, seed)

% [pr_sa_eps, sa_eps_err] = sa_eps(pr, sigma, seed);
%
% Simulate the epsilon contribution to SA to the pseudo-range (PR) 
% and accumulated carrier phase (CPH) measurement error using the RTCA standard
% model of a random bias with zero mean and a standard deviation
% of 23 meters applied to each satellite for the length of the test.
%
% Input:
%   pr          - satellite number and associated PR measurement to corrupt 
%                  with sa_eps errors (m) (nx2) [prn pr]
%   sigma       - standard deviation of the bias to apply (1x1) (optional) 
%                  Default is the RTCA proposed value of 23 meters. 
%   seed        - seed value for random number generator (optional)
%                  Default value is 0. 
%   Note: This sa_eps model does not reflect the slowly varying effect 
%         of the sa_eps error which has periods on the order of a few hours.
% Output:
%   pr_sa_eps  - sa_eps corrupted PR measurements (m) (nx1) 
%   sa_eps_err - sa_eps bias added to PR measurements (m) (nx1)
%
% See also PROPGEPH, SA_CLOCK, CLOCKERR, TROPDLAY

% Written by: Jimmy LaMance 1/15/97
% Copyright (c) 1998 by Constell, Inc.

% Reference: Global Positioning System: Theory and Applications
% Volume 1, Parkinson and Spilker, pages 605.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
p_handle=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on SA_EPS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in optional variables if required.
if nargin < 2
  sigma = 23;        % set to the RTCA default value
end % if nargin > 1
if nargin < 3             % eps_model_data is provided
  seed = 0;     % defualt seed value
end % if nargin == 3 

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'SA_EPS';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = pr;
  
  estruct.variable(2).name = 'sigma';
  estruct.variable(2).req_dim = [1 1];
  estruct.variable(2).var = sigma;
  
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

% allocate the output data arrays
pr_sa_eps = ones(size(pr,1),1) * inf;
sa_eps_err = ones(size(pr,1),1) * inf;

% determine the total number of satellites
% start by sorting the prn numbers and looking for changes
prn_sort = sort(pr(:,1));
prn_change = find(diff(prn_sort) ~= 0);

% create a matrix that has sorted and reduced prns [1 2 4 5 6 8 ... 28]
prn = [prn_sort(prn_change); prn_sort(length(prn_sort))];

% compute the total number of prns 
num_prns = length(prn);

% seed the random number generator
randn('seed',seed);

% generate random number and multiple by the sigma
sv_sa_eps_err = randn(num_prns,1) * sigma;

% loop over the number of satellites and appliy the apporpriate bias to 
% each satellite PR measurement
for i = 1:num_prns
  I = find(pr(:,1) == prn(i));
  if any(I)
    pr_sa_eps(I) = pr(I,2) + sv_sa_eps_err(i);
    sa_eps_err(I) = sv_sa_eps_err(i) * ones(length(I),1);
  else
    fprintf('Warning in processing PR measurements in SA_EPS.\n')
    fprintf('A satellite has been lost in the processing.\n');
    fprintf('Output may be in error.\n\n');
  end % if any(I)
end % for i = 1:num_prns

%%%%% END ALGORITHM CODE %%%%%

% end of SA_EPS
