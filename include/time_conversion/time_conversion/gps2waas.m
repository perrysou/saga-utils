% Function waas_time = gps2waas(gps_time)
% converts gps time, to WAAS system time, which is measured in seconds
% since 1/1/00 12:00:00.
%
% Seebany Datta-Barua
% 11 Mar 2005

function waas_time = gps2waas(gps_time)

waas_time = sec2waas(gpst2sec(gps_time));


