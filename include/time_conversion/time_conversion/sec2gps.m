% Function GPS_TIME = sec2gps(sec) converts
% number of seconds into an n X 2 matrix
% GPS_TIME = [GPS_week_number GPS_time_of_week]
%
% Seebany Datta-Barua
% 18.09.01

function GPS_TIME = sec2gps(sec)

% Define constants
sec_per_day = 60*60*24;
sec_per_wk = sec_per_day * 7;

GPS_week = floor(sec./sec_per_wk);

GPS_tow = mod(sec,sec_per_wk);

GPS_TIME = [GPS_week GPS_tow];

return
