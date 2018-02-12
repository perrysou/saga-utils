% function julian = julian_calculation(t_input)
% This function compute the julian day and julian century from the local
% time and timezone information. Ephemeris are calculated with a delta_t=0
% seconds. 
% Input parameters:
%   time: a structure that specify the time when the sun position is
%   calculated. 
%       time.year: year. Valid for [-2000, 6000]
%       time.month: month [1-12]
%       time.day: calendar day [1-31]
%       time.hour: local hour [0-23]
%       time.min: minute [0-59]
%       time.sec: second [0-59]
%       time.UTC: offset hour from UTC. Local time = Greenwich time + time.UTC
%   This input can also be passed using the Matlab time format ('dd-mmm-yyyy HH:MM:SS'). 
%   In that case, the time has to be specified as UTC time (time.UTC = 0)
%
% This function was cut and pasted from sun_position_sdb.m to make it
% accessible to other functions.
% 
% Seebany Datta-Barua
% 6 June 2008
function julian = julian_calculation(t_input)

% If time input is a Matlab time string, extract the information from
% this string and create the structure as defined in the main header of
% this script.
if ~isstruct(t_input)
    tt = datevec(t_input);
    time.UTC=0;
    time.year = tt(1);
    time.month = tt(2);
    time.day = tt(3);
    time.hour = tt(4);
    time.min = tt(5);
    time.sec = tt(6);
else
    time = t_input;
end

if(time.month == 1 | time.month == 2)
    Y = time.year - 1;
    M = time.month + 12;
else
    Y = time.year;
    M = time.month; 
end
ut_time = ((time.hour - time.UTC)/24) + (time.min/(60*24)) + (time.sec/(60*60*24)); % time of day in UT time. 
D = time.day + ut_time; % Day of month in decimal time, ex. 2sd day of month at 12:30:30UT, D=2.521180556


% In 1582, the gregorian calendar was adopted
if(time.year == 1582)
    if(time.month == 10)
        if(time.day <= 4) % The Julian calendar ended on October 4, 1582
            B = 0;    
        elseif(time.day >= 15) % The Gregorian calendar started on October 15, 1582
            A = floor(Y/100);
            B = 2 - A + floor(A/4);    
        else
            disp('This date never existed!. Date automatically set to October 4, 1582');
            time.month = 10;
            time.day = 4; 
            B = 0;
        end
    elseif(time.month<10) % Julian calendar 
        B = 0;
    else % Gregorian calendar
        A = floor(Y/100);
        B = 2 - A + floor(A/4);
    end
    
elseif(time.year<1582) % Julian calendar
    B = 0;
else
    A = floor(Y/100); % Gregorian calendar
    B = 2 - A + floor(A/4);
end

julian.day = floor(365.25*(Y+4716)) + floor(30.6001*(M+1)) + D + B - 1524.5;

delta_t = 0; % 33.184;
julian.ephemeris_day = julian.day + (delta_t/86400);

julian.century = (julian.day - 2451545) / 36525; 

julian.ephemeris_century = (julian.ephemeris_day - 2451545) / 36525;

julian.ephemeris_millenium = julian.ephemeris_century / 10; 


