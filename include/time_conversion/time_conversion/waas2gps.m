% Function gps_time = waas2gps(waastime)
% converts waas time, measured in seconds
% since 1/1/00 12:00:00 to UTC time.
%
% Seebany Datta-Barua
% 20 Oct 2003
% 23 June 2006 SDB modified to have output dimension nx2 instead of 1x2n.
function gps_time = waas2gps(waastime)

gps_time = sec2gps(waas2sec(waastime));

if size(gps_time,2) ~= 2
    gps_time = reshape(gps_time, length(gps_time)/2, 2);
end