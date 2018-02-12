% Function waastime = sec2waas(sec)
% converts seconds since 1/1/80 start of GPS
% to waas time = seconds since 1/1/00 12:00:00 UT.
%
% Seebany Datta-Barua
% 20 Oct 2003

function waastime = sec2waas(sec)

waastime = sec - 630763213;




