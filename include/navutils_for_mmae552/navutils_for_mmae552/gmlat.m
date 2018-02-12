% Function gmlat, gglat, gglon
% returns geomagnetic latitude given
% geographic latitude and longitude.
% Taken from idlfile gmlat.pro.
%
% Seebany Datta-Barua
% 18 May 2005

function maglat = gmlat(gglat, gglon)

% Magnetic pole location in geographic coordinates
dtor = pi/180;
radeg = 180/pi;

polelat = 78.7*dtor;% deg N converted to rad
polelon = 290.1*dtor;% deg E converted to rad

gglat = gglat*dtor;
gglon = gglon*dtor;

maglat = asin( sin(gglat).*sin(polelat) + ...
	cos(gglat).*cos(polelat).*cos(gglon - polelon) );

maglat = maglat*radeg;