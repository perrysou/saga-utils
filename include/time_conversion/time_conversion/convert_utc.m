function [GPS, ABS, RINEX] = convert_utc(UTC, leap)

% Given the data structure UTC, convert it into
% the other 3 types of time data, GPS, Absolute, Rinex File type
% where the structures contain the following fields :
%
% GPS : week, sec
% ABS : sec
% UTC : year, mon, day, hour, min, sec
% RINEX : year, day
% If leap = 1, then we will account for leap seconds (default)

if (nargin < 2)
    leap = 1;
end;

week2sec = 60 * 60 * 24 * 7;
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31]; % we start in a leap year (1980)

RINEX.year = UTC.year;

curryear = UTC.year;
currmon = UTC.mon;
currday = UTC.day;
totalday = 0;
dayinyear = 0;

if (isleapyear(curryear))
    daysinmonth(2) = 29;
end;

while ((curryear ~= 1980) | ...
       (currmon ~= 1) | ...
       (currday ~= 6))
   totalday = totalday + 1;
   currday = currday - 1;
   % if we're in the same year, keep track of the day of year
   if (curryear == UTC.year)
       dayinyear = dayinyear + 1;
   end;
   
   %keep track of month/year changes
   if (currday == 0)
       currmon = currmon - 1;
       if (currmon == 0)
           curryear = curryear - 1;
           if (isleapyear(curryear))
               daysinmonth(2) = 29;
           else
               daysinmonth(2) = 28;
           end;
           currmon = 12;
       end;
       currday = daysinmonth(currmon);
   end;
end;

if (UTC.year == 1980)
    dayinyear = dayinyear + 6;
end;

RINEX.day = dayinyear;

if (leap)
    delta = leap_correction_utc(UTC);
else 
    delta = 0;
end;
ABS.sec = totalday * 24 * 60 * 60 + UTC.hour * 3600 + UTC.min * 60 + UTC.sec + delta;
GPS.week = floor(ABS.sec / week2sec);
GPS.sec = ABS.sec - GPS.week * week2sec;