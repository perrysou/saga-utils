% Function enu = xyz2enu(xyz,lat,lon,ht)
% takes in an n x 3 matrix xyz for which the columns are
% WGS-84 x,y, and z position, respectively, and
% returns an n x 3 matrix enu of these positions rotated
% into an ENU frame whose origin is at lat, lon (degrees).
%
% Note: requires rot.m and wgslla2xyz.m to be defined
% by adding directories containing these files to the
% Matlab path.
%
% Seebany Datta-Barua
% 20.07.01
% 28.08.01

function enu = xyz2enu(xyz,lat,lon,ht)

refxyz = wgslla2xyz(lat,lon,ht);

R1 = rot(90+lon, 3);
R2 = rot(90-lat, 1);
R = R2*R1;

diffxyz = xyz' - refxyz'*ones([1,size(xyz,1)]);

% ned = ecef2ned(diffxyz', [lat*pi/180, lon*pi/180, ht]);
enu = quicktimes(R,diffxyz);

%for i = 1:size(xyz,1)

%    diffxyz = xyz(i,:)' - refxyz;
%    enu(:,i) = R*diffxyz;
%    enu(:,i) = wgsxyz2enu([xyz(i,1);xyz(i,2);xyz(i,3)],lat,lon,ht,refxyz);
%end

enu = enu';
return

% ----------------------------------------
% Support function
% ----------------------------------------

function Vout = quicktimes(T,Vin)

[m n] = size(T);
if m ~= 3 | n ~= 3
    error('xyz2enu: Rotation matrix must be 3 x 3');
end

[m n] = size(Vin);
if m ~= 3
    error('xyz2enu: Input *vector* must be 3 x n');
end

Vout(1,:) = T(1,1) * Vin(1,:) + T(1,2) * Vin(2,:) + T(1,3) * Vin(3,:);
Vout(2,:) = T(2,1) * Vin(1,:) + T(2,2) * Vin(2,:) + T(2,3) * Vin(3,:);
Vout(3,:) = T(3,1) * Vin(1,:) + T(3,2) * Vin(2,:) + T(3,3) * Vin(3,:);

return
