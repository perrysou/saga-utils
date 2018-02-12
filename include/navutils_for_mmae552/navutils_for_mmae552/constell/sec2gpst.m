function [GPS_time] = sec2gpst(total_gps_secs)

% [GPS_time] = sec2gpst(total_gps_secs);
%
% Utility function to convert gps time in seconds from Jan. 6, 1980 to a
% 2-element matrix of gps_week and gps_sec from beginning of week
%
% Input:  total_gps_seconds - gps_weeks*86400*7 + gps_secs (nx1) (seconds)
%             
% Output: GPS_time - [gps_weeks gps_secs] (nx2)
% 
% Where gps_secs is from the beginning of the week (midnight Saturday night).
%
% See also GPST2SEC

% Written by:  Maria J. Evans   9/7/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
p_handle=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on SEC2GPST for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'SEC2GPST';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'total_gps_secs';
  estruct.variable(1).req_dim = [901 1];
  estruct.variable(1).var = total_gps_secs;
  
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

GPS_time(:,1) = floor(total_gps_secs./(86400*7));
GPS_time(:,2) = total_gps_secs - GPS_time(:,1).*86400*7;

%%%%% END ALGORITHM CODE %%%%%

% end of SEC2GPST
