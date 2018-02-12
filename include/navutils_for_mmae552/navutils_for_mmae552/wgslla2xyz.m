%Function xyz = wgslla2xyz(lat, lon, alt) returns the
%equivalent WGS84 XYZ coordinates (in meters) for a
%given geodetic latitude "lat" (degrees), longitude "lon"
%(degrees), and altitude above the WGS84 ellipsoid
%in meters.  Note: N latitude is positive, S latitude
%is negative, E longitude is positive, W longitude is
%negative.
%
%Ref: Decker, B. L., World Geodetic System 1984,
%Defense Mapping Agency Aerospace Center.
% 30 Sept 2003 Modified for use instead of llh2xyz.m SDB.
% 1 Feb 2005 Brought into service on Windows machine because no copy of
% llh2xyz.m has been ported over here.
% 11 Nov 2005 Adjusting out of range lats and lons.
% 30 May 2013 Altitude in meters.

function xyz = wgslla2xyz(wlat, wlon, walt)


rows = find(wlon > 180);
wlon(rows) = wlon(rows) - 360;
rows = find(wlon < -180);
wlon(rows) = wlon(rows)+360;
rows = find(wlat > 90);
if ~isempty(rows)
    keyboard
end
%                 llh2(rows,2) = llh2(rows,2)+360;


A_EARTH = 6378137;
flattening = 1/298.257223563; 
NAV_E2 = (2-flattening)*flattening; % e^2
deg2rad = pi/180;

slat = sin(wlat.*deg2rad);
clat = cos(wlat.*deg2rad);
r_n = A_EARTH./sqrt(1 - NAV_E2.*slat.*slat);
xyz = zeros(length(wlat),3);
xyz(:,1) = (r_n + walt).*clat.*cos(wlon*deg2rad);
xyz(:,2) = (r_n + walt).*clat.*sin(wlon*deg2rad);
xyz(:,3) = (r_n.*(1 - NAV_E2) + walt).*slat;

%size(xyz)

if ((wlat < -90.0) | (wlat > +90.0) |...
        (wlon < -180.0) | (wlon > +360.0))
    error('WGS lat or WGS lon out of range');
end
return


