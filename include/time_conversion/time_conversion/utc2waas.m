% Function waastime = utc2waas(utc_time)
% converts UTC time = [yy mm dd hh min ss]
% to waas time = seconds since 1/1/00 12:00:00 UT.
%
% Seebany Datta-Barua
% 20 Oct 2003

function waastime = utc2waas(utc_time)

waastime = sec2waas(gpst2sec(utc2gps(utc_time)));



