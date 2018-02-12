function [total_gps_seconds] = gpst2sec(GPS_time)

% [total_gps_seconds] = gpst2sec(GPS_time)
%
% Utility function to convert a GPS time structure to pure seconds.  This is
% useful when generating time intervals that may cross a week boundary.
%
%   Input:  GPS_time = [gps_weeks gps_secs (from beginning of week)]
%            valid GPS_week values are 1-3640 (years 1980-2050)
%            valid GPS_sec values are 0-604799
%
%   Output: total_gps_seconds = gps_weeks*86400*7 + gps_secs
%
% See also SEC2GPST

% Written by:  Maria J. Evans  9/7/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
UTC_time=[]; leap_sec=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on GPST2SEC for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'GPST2SEC';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'GPS_time';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = GPS_time;
  estruct.variable(1).type = 'GPS_TIME';
  
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

total_gps_seconds = GPS_time(:,1).*86400*7 + GPS_time(:,2);

%%%%% END ALGORITHM CODE %%%%%

% end of GPST2SEC


