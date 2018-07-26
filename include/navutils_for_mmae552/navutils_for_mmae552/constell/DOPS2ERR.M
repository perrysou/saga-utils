function [pos_err] = dops2err(dops,sigma_pr);

% [pos_err] = dops2err(dops,sigma_pr);
%
% Converts from Dilution of Precision (DOP) to equilivant position error.
% This is a first order approximation that takes the DOP value and multiplies
% by the sigma on the PR measurement to provide an estimate of the position
% error components in the same reference frame as the DOP values.
%
% Input:
%   dops     - DOP values (nx5) [GDOP PDOP HDOP VDOP TDOP]
%   sigma_pr - sigma on the pseudo-range measurements (1x1) (m) (optional)
%               default = 30. 
% Output:
%   pos_err  - estaimte of the position error based on the DOPs and the
%               sigma on the pseudo-range (nx3) (m)
%               [tot_pos_err horizontal_err vertical_err]
%
% See also NED2DOPS, LL2DOPS, PSEUDO_R, DIFFCORR

% Written by: Jimmy LaMance 9/3/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
pos_err=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on DOPS2ERR for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Fill in the optional variables if not included in the input arguments
if nargin < 2
  sigma_pr = 30;
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'DOPS2ERR';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'dops';
  estruct.variable(1).req_dim = [901 5];
  estruct.variable(1).var = dops;
  
  estruct.variable(2).name = 'sigma_pr';
  estruct.variable(2).req_dim = [1 1];
  estruct.variable(2).var = sigma_pr;
  
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

% allocate the output data
pos_err = ones(size(dops,1),3) * inf;

% fill in the pos_err matrix
pos_err(:,1) = dops(:,2) * sigma_pr;
pos_err(:,2) = dops(:,3) * sigma_pr;
pos_err(:,3) = dops(:,4) * sigma_pr;

%%%%% END ALGORITHM CODE %%%%%

% end of DOPS2ERR
