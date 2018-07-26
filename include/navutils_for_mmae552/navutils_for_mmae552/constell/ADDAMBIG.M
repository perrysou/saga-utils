function [cphnew, ambig] = addambig(cph, seed);

% [cphnew, ambig] = addambig(cph, seed);
%
% Adds a random set of integer ambiguities (N) to the accumulated 
% carrier phase (CPH) measurements.  The ambiguities will between -1e6 and 1e6.
%
% Input:
%   cph     - satellite number and associated CPH measurement to apply
%              N ambiguities (cycles) (nx2) [prn cph]
%   seed    - seed value for random number generator (optional)
%              Default value is 0. 
% Output:
%   cphnew  - CPH measurements with N integer ambiguities added (cycles) (nx1) 
%   ambig   - N integer ambiguities (nx1)
%
% See also PSEUDO_R, SA_CLOCK, CLOCKERR, TROPDLAY

% Written by: Jimmy LaMance 9/15/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
cphnew=[]; ambig=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ADDAMBIG for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check if a seed is provided
if nargin < 2
  seed = 0;
end % if nargin < 6

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ADDAMBIG';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'cph';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = cph;
  
  estruct.variable(2).name = 'seed';
  estruct.variable(2).req_dim = [1 1];
  estruct.variable(2).var = seed;
  
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
cphnew = ones(size(cph,1),1) * inf;
ambig = ones(size(cph,1),1) * inf;

% determine the total number of satellites
% start by sorting the prn numbers and looking for changes
prn_sort = sort(cph(:,1));
prn_change = find(diff(prn_sort) ~= 0);

% create a matrix that has sorted and reduced prns [1 2 4 5 6 8 ... 28]
prn = [prn_sort(prn_change); prn_sort(length(prn_sort))];

% compute the total number of prns 
num_prns = length(prn);

% seed the random number generator
randn('seed',seed);

% generate random number and multiple by the sigma
cph_amb = randn(num_prns,1) * 1e8;

% change the sign on a random set of these ambiguities
sign_chng = randn(num_prns,1);

I_sign_change = find(sign_chng > .5);

if ~isempty(I_sign_change) 
  cph_amb(I_sign_change) = -cph_amb(I_sign_change);
end % if ~isempty(I_sign_change) 

% loop over the number of satellites and appliy the apporpriate bias to 
% each satellite PR measurement
for i = 1:num_prns
  I = find(cph(:,1) == prn(i));
  if any(I)
    cphnew(I) = cph(I,2) + cph_amb(i);
    ambig(I) = ones(length(I),1) * cph_amb(i);
  else
    fprintf('Warning in processing CPH measurements in ADDAMBIG.\n')
    fprintf('A satellite has been lost in the processing.\n');
    fprintf('Output may be in error.\n\n');
  end % if any(I)
end % for i = 1:num_prns

%%%%% END ALGORITHM CODE %%%%%

% end of ADDAMBIG
