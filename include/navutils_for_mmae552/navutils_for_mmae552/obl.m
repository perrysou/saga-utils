% Function obliquity = obl(el, alt_H_iono) returns the obliquity factor in the range
% of [1,5] to the caller.  Assumes thin shell height of ionosphere is 350
% km if 2nd input alt_H_iono is not specified, in meters.  Input argument 
% el should be in radians.
% 
% Seebany Datta-Barua
% 12 Nov 2003
% 4 June 2005 Modified to allow variable shell height.

function obliquity = obl(el, alt_H_iono)

if nargin == 1
    obliquity = 1 ./ cos(asin(0.94797966*cos(el)));
else
    gpsconst;
    obliquity = 1 ./ sqrt( 1 - (R_e*cos(el)./(R_e + alt_H_iono) ).^2 ) ;
end