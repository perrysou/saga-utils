% function [gclat, gcrad] = gd2gc(gdlat, gdht)
% converts geodetic (i.e. WGS 84) latitude to geocentric (i.e. spherical
% earth) latitude. Latitudes gdlat and gclat in degrees, height gdht and
% spherical radius gcrad in meters.
%
% Seebany Datta-Barua
% 14 Sept 2006
function [gclat, gcrad] = gd2gc(gdlat, gdht)

load gpsconst R_e

gdlatrad = gdlat*pi/180;

% Flattening taken from wgslla2xyz.m.
f = 1/298.257223563;
NAV_E2 = (2-f)*f; % also e^2
r_n = R_e./sqrt(1 - NAV_E2.*(sin(gdlatrad)).^2);


% Geocentric latitude at surface.
lambda_s = atan( (1-f)^2*tan(gdlatrad));

% Geocentric latitude away from surface.
gclat = atan2( (gdht*sin(gdlatrad) + R_e*sin(lambda_s)), ...
    (gdht*cos(gdlatrad) + R_e*cos(lambda_s)) );

% Expression for spherical radius taken from Bate, Mueller, White, p. 95.
x = (r_n + gdht).*cos(gdlatrad);
z = (r_n*(1-NAV_E2) + gdht).*sin(gdlatrad);
gcrad = sqrt(x.^2 + z.^2);

gclat = gclat*180/pi;
