function [GPS, UTC, RINEX] = convert_abs(ABS, leap)

% Given the data structure ABS (aboslute time), convert it into
% the other 3 types of time data, GPS, UTC, Rinex File type
% where the structures contain the following fields :
%
% GPS : week, sec
% ABS : sec
% UTC : year, mon, day, hour, min, sec
% RINEX : year, day
% If leap = 1, then we will account for leap seconds (default)
%
% 29 June 2006 Discovered bug: It only returns one date, even if you pass a
% vector of times spanning multiple days. SDB

if (nargin < 2)
    leap = 1;
end;

week2sec = 60 * 60 * 24 * 7;
daysinmonth = [31 29 31 30 31 30 31 31 30 31 30 31]; % we start in a leap year (1980)

GPS.week = floor(ABS.sec / week2sec);
GPS.sec = ABS.sec - GPS.week * week2sec;

if (leap)
    delta = leap_correction_gps(ABS.sec);
    secondtotal = ABS.sec - delta;
else
    secondtotal = ABS.sec;
end;

totaldays = floor(secondtotal / (24 * 60 * 60));
leftover = secondtotal - totaldays * 24 * 60 * 60;
UTC.hour = floor (leftover / 3600);
leftover = leftover - UTC.hour * 3600;
UTC.min = floor (leftover / 60);
UTC.sec = leftover - UTC.min * 60;

curryear = 1980;
currmon = 1;
currday = 6;
dayinyear = 6;

for i = 1:totaldays
    currday = currday + 1;
    dayinyear = dayinyear + 1;
    if (currday == daysinmonth(currmon) + 1)
        currday = 1;
        currmon = currmon + 1;
    end;
    if (currmon == 13)
        currmon = 1;
        dayinyear = 1;
        curryear = curryear + 1;
        if (isleapyear(curryear))
            daysinmonth(2) = 29;
        else
            daysinmonth(2) = 28;
        end;
    end;
end;

UTC.year = curryear;
UTC.mon = currmon;
UTC.day = currday;
RINEX.day = dayinyear;
RINEX.year = curryear;