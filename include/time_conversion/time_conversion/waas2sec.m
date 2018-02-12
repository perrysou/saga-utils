% Function sec = waas2sec(waastime)
% converts waas time, measured in seconds
% since 1/1/00 12:00:00 to total seconds
% since 1/1/80 when GPS time started.
%
% Seebany Datta-Barua
% 20 Oct 2003

function sec = waas2sec(waastime)

sec = waastime + 630763213;
