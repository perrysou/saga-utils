% Function enu = wgsxyz2enu(xyz, lat,lon,ht,refxyz) returns
% a 3 x 1 vector enu which represents the East, North, and Up
% coordinates (in meters) of a point with WGS84 xyz coordinates
% xyz (3 x 1 vector, units in meters) in an ENU coordinate system
% located at wgsxyz coords refxyz (equivalent to wgs coords lat, lon, ht).
%
% Note: Requires function rot.m to be in the 
% same directory

function enu=wgsxyz2enu(xyz,lat,lon,ht,refxyz)

[m n] = size(xyz);
if m ~= 3 | n ~= 1
	error('wgsxyz2enu: xyz input vector must be 3 x 1');
end

[m n] = size(lat);
if m ~= 1 | n ~= 1
	error('wgsxyz2enu: lat input vector must be scalar');
end

[m n] = size(lon);
if m ~= 1 | n ~= 1
	error('wgsxyz2enu: lon input vector must be scalar');
end

[m n] = size(ht);
if m ~= 1 | n ~= 1
	error('wgsxyz2enu: ht input vector must be scalar');
end

[m n] = size(refxyz);
if m ~= 3 | n ~= 1
    error('wgsxyz2enu: refxyz input vector must be 3 x 1');
end


% Difference xyz from reference point
diffxyz = xyz - refxyz;

% Now rotate the (often short) diffxyz vector to enu frame

R1=rot(90+lon, 3);
R2=rot(90-lat, 1);
R=R2*R1;

enu=R*diffxyz;

return;

           
