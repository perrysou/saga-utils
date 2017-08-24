% function [ipplat, ipplon] = sill(az, el, stnlat, stnlon, stnht, alt_H_iono) returns the
% sub-ionospheric latitude and longitude, in degrees, of the IPP assuming 
% a fixed shell height, given vectors 0<= az<= 2pi and el in radians,
% vectors stnlat and stnlon in degrees, and stnht and alt_H_iono in m.
%
% Adapted from the sill.pro IDL function supplied by Pat Doherty.
% Seebany Datta-Barua
% 5 May 2004
% 10 Aug 2006 Making alt_H_iono mandatory for use with 1300 km iono end ht.
% 23 Aug 2006 Adding input stnht for use with space-based GPS.
% 3 Sep 2006 Issues with dlon being imag because of rounding error acos > 1.
% 8 Dec 2006 acos > 1 means there's no solution, by law of sines.  Return
% nans.
% 22 Oct 2007 Testing Hiro's derivation by rotations.


function [ipplat, ipplon] = sill(az, el, stnlat, stnlon, stnht, alt_H_iono)
gpsconst;
% if nargin == 5
    H_iono = alt_H_iono;
% end

ipplat = zeros(length(az),1);
ipplon = zeros(length(az),1);

deg2rad = pi/180;
rad2deg = 180/pi;
if size(stnlat,1) == 1
    stnlat = ones(size(az))*stnlat;
    stnlon = ones(size(az))*stnlon;
    stnht = ones(size(az))*stnht;
end

stnlat = deg2rad*stnlat;
stnlon = deg2rad*stnlon;

% Radius of earth at station
Re = R_e/1e3 * ( .99832 + .00168*cos(2*stnlat));

% In GPS textbook, chi = xi'.
% If (Re + stnht/1e3).* cos(el) > (Re+H_iono/1e3) there is no triangle for
% the law of sines to work.  Replace with nans.
chi = asin( cos(el) .* (Re + stnht/1e3) ./ (Re + H_iono/1e3) );
rows = find((Re + stnht/1e3).* cos(el) > (Re+H_iono/1e3));
chi(rows) = nan;

% A is earth-centered angle in rad
A = pi/2 - chi - el;
ipplat = asin( sin(stnlat) .*cos(A) + cos(stnlat) .* sin(A) .* cos(az));

% Rounding errors make acos imaginary.  Force to be real.
% Changed the rounding error to tolerate a higher value than 1e-15. SDB
% 130903.
if max(( cos(A) - sin(ipplat).*sin(stnlat) ) ./ ( cos(ipplat).*cos(stnlat) )) -1>1e-14
    ipplat = nan(size(az));
    ipplon = nan(size(az));
    warning('sill.m line 48: acos > 1 means no solution.  Returning nans.')
    return
else
dlon = acos ( ( cos(A) - sin(ipplat).*sin(stnlat) ) ./ ( cos(ipplat).*cos(stnlat) ) );
    rows = find(imag(dlon));
    dlon(rows) = 0;
end
clear rows
rows = find(az <= pi);
if ~isempty(rows)
% rows = find(dlon);
    ipplon(rows) = stnlon(rows) + dlon(rows);
end

rows = find(az > pi);
if ~isempty(rows)
%     keyboard
    ipplon(rows) = stnlon(rows) - dlon(rows);
end

ipplat = rad2deg*ipplat;
ipplon = rad2deg*ipplon;

% % Alternate method of getting ipplon:
% R1 = [1 0 0; 0 cos(A) sin(A); 0 -sin(A) cos(A)];
% R2 = [cos(az) sin(az) 0; -sin(az) cos(az) 0; 0 0 1];
% R3 = [1 0 0; 0 cos(stnlat) sin(stnlat); 0 -sin(stnlat) cos(stnlat)];
% R = R3*R2*R1;
% xyz = R*[0;0;1];
% ipplat2 = asin(xyz(2))*180/pi
% dlon2 = atan2(xyz(1),xyz(3));
% ipplon2 = (dlon2+stnlon)*180/pi
% 
% % Test the rotation back to xyz frame.
% R4 = [cos(-stnlon) 0 -sin(-stnlon); 0 1 0; sin(-stnlon) 0 cos(-stnlon)];
% xyz = R4*R*[0;0;1];
% ipplat3 = asin(xyz(2))*180/pi
% ipplon3 = atan2(xyz(1),xyz(3))*180/pi