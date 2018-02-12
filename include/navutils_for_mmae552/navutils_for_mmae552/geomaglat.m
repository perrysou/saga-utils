% function gmlat = geomaglat(gclat,gclon)
% returns the geomagnetic latitude of a point
% with geocentric (geographic) coords gclat,gclon.
% Note that gclon must be in positive degrees 
% east of 0.  Based on values given by Sam Pullen 
% and available at 
% http://nssdc.gsfc.nasa.gov/space/cgm/cgm.html
%
% Seebany Datta-Barua
% 26.06.01

function gmlat = geomaglat(gclat,gclon)

if (gclon >= 180)
    gclon = gclon-360;
elseif (gclon < -180)
    gclon = gclon + 360;
end

if (gclon > 20 & gclon < 180)
    error('Coordinates given are beyond range of available geomagnetic data');
end

% Load table of location of gm equator in gccoords, 
% as a function gclon.  This is the data given by Sam.  
% First column is longitude.  Second column is latitude 
% of geomagnetic equator relative to geographic equator 
% at that longitude.
load gmeq.mat

if mod(gclon,10) == 0
    gmlat = gclat - gmeq(find(gmeq(:,1) == gclon),2);
else
    % Linearly interpolate offsets in gmlat for 
    % longitudes between the values given in 
    % Sam's table.
    x1 = 10*floor(gclon/10);
    x2 = 10*ceil(gclon/10);
    y1 = gmeq(find(gmeq(:,1) == x1),2);
    y2 = gmeq(find(gmeq(:,1) == x2),2);

    % Find equation of the line y = mx+b.
    m = (y2-y1)/(x2-x1);
    b = y2 - m*x2;

    gmlat = gclat - (m*gclon + b);
end
return









