function [gst] = sidereal(UTC_time)

% [gst] = sidereal(UTC_time);
% 
% This function computes the Greenwich sidereal time (GST).
% GST is defined as the angle (or time converted to angle) between the
% inertial x-axis, xi (first point of Aries), and the Greenwich meridian.
%
% Input:  
%   UTC_time - matrix of the form [year month day hour minute second] (nx6)
%               with 4-digit (1980) years, valid years are 1900 - 2079 
% Output: 
%   gst     - angle between the ECI and ECEF coordinate system (rad) (nx1)
%
% See also ECI2ECEF, ECEF2ECI, GPS2UTC, UTC2GPS

% Written by: Jimmy LaMance 12/9/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, UTC2GPS

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
gst=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on SIDEREAL for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'SIDEREAL';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'UTC_time';
  estruct.variable(1).req_dim = [901 6];
  estruct.variable(1).var = UTC_time;
  
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

jd_gps = 2444244.5; % Julian date on Jan. 6, 1980, 0:0:0 GMT 
jd0 = 2415020.;     % Julian date on Jan. 0, 1900, 12 noon GMT 

% find the number of days since the GPS epoch
% call utc2gps so we can compute the number of GPS days using existing functions
[GPS_week GPS_sec GPS_day] = utc2gps(UTC_time,0);

% compute the julian date (JD_GPS_EPOCH + days since JD0)
julian_date = jd_gps + fix(GPS_day);
julian_century =  (julian_date - jd0) ./ 36525.;    

gst = 99.6909833 + 36000.7689 .* julian_century + ...
       0.00038708 .* julian_century.^2 + (GPS_day-fix(GPS_day))*1440.*.25068447;

gst = rem(gst,360.0);      % want gst: 0 -> 360 deg

gst = gst * pi / 180;      % return gst in radians

%%%%% END ALGORITHM CODE %%%%%

% end of SIDEREAL 
