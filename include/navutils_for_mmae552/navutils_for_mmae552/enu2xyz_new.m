% Function xyz = enu2xyz_new(enu, lat, lon, ht) takes in an n x 3 matrix enu
% for which the cols are in an ENU frame relative to the origin given by
% lat, lon, ht (degrees, m).  The matrix returned xyz, is the WGS-84 x,y,z
% position.
%
% Note: requires rot.m and wgslla2xyz.m to be defined.
%
% Seebany Datta-Barua
% 14 Nov 2003
% 8 Aug 2006 Modifying to accept vector lat, lon, ht.

function xyz = enu2xyz_new(enu, lat, lon, ht)

refxyz = wgslla2xyz(lat,lon,ht);
% refxyz = llh2xyz([lat,lon,ht]);

% Rotation from xyz to enu.  Need to invert this.
R1=rot_new(90+lon, 3);
R2=rot_new(90-lat, 1);

if length(lat) == 1
R = R2*R1;
diffxyz = quicktimes(inv(R),enu');
% 15 July 2004.  SDB. Had to add "'" to refxyz to get dimensions to work on
% iono_map_malicious_km applied to 031122 data.
xyz = diffxyz + refxyz'*ones([1,size(enu,1)]);
else
    for i = 1:length(lat)
        R = R2(:,:,i)*R1(:,:,i);        
        diffxyz(:,i) = inv(R)*enu(i,:)';
    end
    xyz = diffxyz' + refxyz;
end