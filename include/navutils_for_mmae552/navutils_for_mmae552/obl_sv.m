% Function obliquity = obl_sv(el, H_sv, H_iono) returns the obliquity factor in the range
% of [1,5] to the caller.  H_sv is altitude of receiver and H_iono is
% altitude of ionosphere shell, in meters.  Input argument
% el should be in radians.
%
% Seebany Datta-Barua
% 12 Nov 2003
% 4 June 2005 Modified to allow variable shell height.
% 14 July 2006 Modified to allow for space-borne receivers.  Don't know if
% this works for elevation < 0.
% 8 Dec 2006 For el < 0, there exist two triangles (by law of sines) if 
% (R_e + H_sv).*cos(el) < (R_e + alt_H_iono), one for equality, and none for >.

function obliquity = obl_sv(el, H_sv, alt_H_iono)

gpsconst;
obliquity = 1 ./ sqrt( 1 - ( ...
    (R_e+H_sv).*cos(el)./(R_e + alt_H_iono) ).^2 ) ;

% No solutions.
[r,c] = find((R_e+H_sv).*cos(el)./(R_e + alt_H_iono) > 1);
obliquity(r,c) = nan;
