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
% 12 May 2006 Modifying to accept vector lat, lon, ht.  SDB

function enu = xyz2enu_new(xyz,lat,lon,ht)

refxyz = wgslla2xyz(lat,lon,ht);

R1 = rot_new(90+lon, 3);
R2 = rot_new(90-lat, 1);

if length(lat) == 1
R = R2*R1;
diffxyz = xyz' - refxyz'*ones([1,size(xyz,1)]);
enu = quicktimes(R,diffxyz);
else
    diffxyz = xyz' - refxyz';
    for i = 1:length(lat)
        R = R2(:,:,i)*R1(:,:,i);
        enu(:,i) = R*diffxyz(:,i);
    end
end
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
