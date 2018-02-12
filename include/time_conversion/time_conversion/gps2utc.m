function [UTC_time, leap_sec, day_of_year] = gps2utc(GPS_time, offset)

% [UTC_time, leap_sec] = gps2utc(GPS_time, offset); 
%
% Converts a GPS time into an equivalent UTC time (matrices). 
%
% Input:  
%   GPS_time - GPS time (nx2) [GPS_week GPS_sec]
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               GPS week values are kept in linear time without the 1024 rollovers
%               (e.g. UTC time [1999 10 1 0 0 0] = GPS time [1029 432013])
%               use the MATLAB "mod" function to obtain 1024 based weeks
%                  new_GPS_week = mod(GPS_week,1024);
%               valid GPS_sec values are 0-604799
%   offset   - leap seconds for the GPS times (optional) (1x1 or nx1)
%               valid offset values are 0-500
% Output: 
%   UTC_time - matrix of the form [year month day hour minute second]
%               with 4-digit year (1980), nx6 matrix
%   leap_sec - leap seconds applied to UTC relative to GPS (optional)
%
% Note: Any invalid inputs times will result in the UTC equilivant time
%       being filled with inf (infinity) and a warning will be issued.  If
%       all of the GPS time input is invalid, the function will terminate
%       with an error.
%
% See also UTC2GPS, GPS2LEAP

% Written by: Maria Evans/Jimmy LaMance 10/16/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, UTC2LEAP

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
UTC_time=[]; leap_sec=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on GPS2UTC for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

if nargin < 2
  offset = zeros(size(GPS_time,1),1);
end % if nargin < 2

if nargin == 2 & size(offset,1) == 1 & size(offset,2) == 1
  offset = offset * ones(size(GPS_time,1),1);
end
 
% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'GPS2UTC';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'GPS_time';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = GPS_time;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'offset';
  estruct.variable(2).req_dim = [901 1; 1 1];
  estruct.variable(2).var = offset;
                                
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

% Break out GPS week and seconds into separate variable
GPS_week = GPS_time(:,1);
GPS_sec = GPS_time(:,2);

% allocate the momory for the UTC_time working matrix
UTC_time = ones(size(GPS_week,1),6) * inf;

% compute gpsday and gps seconds since start of GPS time 
gpsday = GPS_week * 7 + GPS_sec ./ 86400;
gpssec = GPS_week * 7 * 86400 + GPS_sec;

% get the integer number of days
total_days = floor(gpsday);

% temp is the number of completed years since the last leap year (0-3)
% the calculation begins by computing the number of full days since
% the beginning of the last leap year.  This is accomplished through
% the rem statement.  Since GPS time started at
% 00:00 on 6 January 1980, five days must be added to the total number
% of days to ensure that the calculation begins at the beginning of a
% leap year.  By subtracting one from this result, the extra day in
% the first year is effectively removed, and the calculation can
% simply be computed by determining the number of times 365 divides
% into the number of days since the last leap year.  On the first day
% of a leap year, the result of this calculation is -1
% so the second statement is used to trap this case.

temp = floor((rem((total_days+5),1461)-1) ./ 365);
I_temp=find(temp < 0);
if any(I_temp), 
  temp(I_temp) = zeros(size(temp(I_temp))); 
end % if

% compute the year
UTC_time(:,1) = 1980 + 4 * floor((total_days + 5) ./ 1461) + temp;

% data matrix with the number of days per month for searching 
% for the month and day
% days in full months for leap year
leapdays =   [0 31 60 91 121 152 182 213 244 274 305 335 366];  
% days in full months for standard year
noleapdays = [0 31 59 90 120 151 181 212 243 273 304 334 365];                                                      

% Leap year flag
% determine which input years are leap years
leap_year = ~rem((UTC_time(:,1)-1980),4);     
I_leap = find(leap_year == 1);                % find leap years
I_no_leap = find(leap_year == 0);             % find standard years

% establish the number of days into the current year
% leap year
if any(I_leap)
  day_of_year(I_leap) = rem((total_days(I_leap) + 5),1461) + 1;                        
end % if any(I_leap)

% standard year
if any(I_no_leap)
  day_of_year(I_no_leap) = ...
      rem(rem((total_days(I_no_leap) + 5),1461) - 366, 365) + 1;  
end % if any(I_no_leap)

% generate the month, loop over the months 1-12 and separate out leap years
for iii = 1:12
  if any(I_leap)
    I_day = find(day_of_year(I_leap) > leapdays(iii));
    UTC_time(I_leap(I_day),2) = ones(size(I_day')) * iii;
    clear I_day
  end % if any(I_leap) 
  
  if any(I_no_leap)
    I_day = find(day_of_year(I_no_leap) > noleapdays(iii));
    UTC_time(I_no_leap(I_day),2) = ones(size(I_day')) * iii;
    clear I_day
  end % if any(I_no_leap)
end % for

% use the month and the matrix with days per month to compute the day 
if any(I_leap)
  UTC_time(I_leap,3) = day_of_year(I_leap)' - leapdays(UTC_time(I_leap,2))';
end % if any(I_leap)

if any(I_no_leap)
  UTC_time(I_no_leap,3) = day_of_year(I_no_leap)' - ...
                          noleapdays(UTC_time(I_no_leap,2))';
end % if any(I_no_leap)

% compute the hours
fracday = rem(GPS_sec, 86400);              % in seconds!

UTC_time(:,4) = fix(fracday ./ 86400 .* 24);

% compute the minutes 
UTC_time(:,5) = fix((fracday - UTC_time(:,4) .* 3600) ./ 60 );

% compute the seconds
UTC_time(:,6) = fracday - UTC_time(:,4) .* 3600 - UTC_time(:,5) .* 60;

% Compensate for leap seconds
% check the input agrument list for offset (leap seconds)
if nargin < 2         % offset is not given and must be computed
  % Call utc2leap to compute the leap second offset for each time
  leap_sec = utc2leap(UTC_time);
else
  leap_sec = offset;
end % if num_args < 2

UTC_time(:,6) = UTC_time(:,6) - leap_sec;

% Check to see if leap_sec offset causes a negative number of seconds
I_shift = find(UTC_time(:,6) < 0);
UTC_time(I_shift,5) = UTC_time(I_shift,5) - 1;
UTC_time(I_shift,6) = UTC_time(I_shift,6) + 60;

% Check to see if the leap second offset causes a negative number of minutes
I_shift = find(UTC_time(:,5) < 0);
UTC_time(I_shift,4) = UTC_time(I_shift,4) - 1;
UTC_time(I_shift,5) = UTC_time(I_shift,5) + 60;

% Check to see if the leap second offset causes a negative number of hours
I_shift = find(UTC_time(:,4) < 0);
UTC_time(I_shift,3) = UTC_time(I_shift,3) - 1;
UTC_time(I_shift,4) = UTC_time(I_shift,4) + 24;

% Check to see if this causes a 0 day value
I_shift = find(UTC_time(:,3) <= 0);
UTC_time(I_shift,2) = UTC_time(I_shift,2) - 1;
I_yr_shift = find(UTC_time(:,2) <= 0);
UTC_time(I_yr_shift,1) = UTC_time(I_yr_shift,1) - 1;
UTC_time(I_yr_shift,2) = UTC_time(I_yr_shift,2) + 12;

% Leap year flag
 % determine which input years are leap years
leap_year = ~rem((UTC_time(I_shift,1)-1980),4);    
I_leap = find(leap_year == 1);                % find leap years
I_no_leap = find(leap_year == 0);             % find standard years

if any(I_leap),
  UTC_time(I_shift(I_leap),3) = leapdays(UTC_time(I_shift(I_leap),2) + 1)' ...
                               -leapdays(UTC_time(I_shift(I_leap),2))';
end;
if any(I_no_leap),
  UTC_time(I_shift(I_no_leap),3) = ...
         noleapdays(UTC_time(I_shift(I_no_leap),2) + 1)' ...
         -noleapdays(UTC_time(I_shift(I_no_leap),2))';
end;

day_of_year = day_of_year';

%%%%% END ALGORITHM CODE %%%%%

% end of GPS2UTC
