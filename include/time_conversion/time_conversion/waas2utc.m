% Function utc_time = waas2utc(waastime)
% converts waas time, measured in seconds
% since 1/1/00 12:00:00 to UTC time.
%
% Seebany Datta-Barua
% 20 Oct 2003

function utc_time = waas2utc(waastime)

utc_time = gps2utc(sec2gps(waas2sec(waastime)));


