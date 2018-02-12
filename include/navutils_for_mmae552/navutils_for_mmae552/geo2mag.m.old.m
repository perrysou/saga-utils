% Function geo2mag.m takes geographic nx1 vectors glat, glon (in degrees) 
% as arguments and returns the geomagnetic lat, lon.  Adapted from 
% geo2mag.pro, downloaded from
% http://astro.uni-tuebingen.de/software/idl/astrolib/astro/geo2mag.html
%
% Seebany Datta-Barua
% 19 May 2005

% NAME:
%       GEO2MAG()
%
% PURPOSE:
%       Convert from geographic to geomagnetic coordinates
% EXPLANATION:
%       Converts from GEOGRAPHIC (latitude,longitude) to GEOMAGNETIC (latitude,
%       longitude).   (Altitude remains the same)
%
%       Latitudes and longitudes are expressed in degrees.
%
% CALLING SEQUENCE:
%       mcoord=geo2mag(gcoord)
%
% INPUT:
%       gcoord = a 2-element array of geographic (latitude,longitude), or an
%                array (2,n) of n such coordinates.
%
% KEYWORD INPUTS:
%       None
%
% OUTPUT:
%       a 2-element array of magnetic (latitude,longitude), or an array (2,n)
%         of n such coordinates
%
% COMMON BLOCKS:
%       None
%
% EXAMPLES:
%       geographic coordinates of magnetic south pole
%
%       IDL> mcoord=geo2mag((79.3,288.59))
%       IDL> print,mcoord
%       89.999992      -173.02325
%
% MODIFICATION HISTORY:
%       Written by Pascal Saint-Hilaire (Saint-Hilaire@astro.phys.ethz.ch),
%            May 2002
%
%-

%====================================================================================
function [mlat, mlon] = geo2mag(glat,glon)

% SOME 'constants'...
Dlon=288.59;   % longitude (in degrees) of Earth's magnetic south pole
%(which is near the geographic north pole!) (1995)
Dlat=79.30;     % latitude (in degrees) of same (1995)
R = 1;          % distance from planet center (value unimportant --
%just need a length for conversion to rectangular coordinates)

% convert first to radians
Dlon=Dlon*pi/180;
Dlat=Dlat*pi/180;

glat=glat*pi/180;
glon=glon*pi/180;
galt=glat * 0. + R; % Creates vector alt, just to get dimensions right...
% in matrix multiplication.

% Make a 3xn matrix.
coord=[glat';glon';galt'];

%convert to rectangular coordinates
%       X-axis: defined by the vector going from Earth's center towards
%            the intersection of the equator and Greenwitch's meridian.
%       Z-axis: axis of the geographic poles
%       Y-axis: defined by Y=Z^X
x=coord(3,:).*cos(coord(1,:)).*cos(coord(2,:));
y=coord(3,:).*cos(coord(1,:)).*sin(coord(2,:));
z=coord(3,:).*sin(coord(1,:));
% xyz = wgslla2xyz(coord(1),coord(2),coord(3))


%Compute 1st rotation matrix : rotation around plane of the equator,
%from the Greenwich meridian to the meridian containing the magnetic
%dipole pole.
glon2mlon=zeros(3,3);
glon2mlon(1,1)=cos(Dlon);
glon2mlon(1,2)=sin(Dlon);
glon2mlon(2,1)=-sin(Dlon);
glon2mlon(2,2)=cos(Dlon);
glon2mlon(3,3)=1;
% out=geolon2maglon * [x,y,z];
out = quicktimes(glon2mlon,[x;y;z]);

%Second rotation : in the plane of the current meridian from geographic
%                  pole to magnetic dipole pole.
tomaglat=zeros(3,3);
tomaglat(1,1)=cos(pi/2-Dlat);
tomaglat(1,3)=-sin(pi/2-Dlat);
tomaglat(3,1)=sin(pi/2-Dlat);
tomaglat(3,3)=cos(pi/2-Dlat);
tomaglat(2,2)=1;
% out= tomaglat * out;
out = quicktimes(tomaglat,out);

%convert back to latitude, longitude and altitude
mlat=atan2(out(3,:),sqrt(out(1,:).^2+out(2,:).^2));
mlat=mlat'*180./pi;
mlon=atan2(out(2,:),out(1,:));
mlon=mlon'*180./pi;
%malt=sqrt(out(:,1)^2+out(:,2)^2+out(:,3)^2)-R
%  I don't care about that one...just put it there for completeness' sake

%===============================================================================
