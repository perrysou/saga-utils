% function mENU = genu2menu(lat, lon, ht, yyyy, doy, gENU)
% rotates a vector gENU = [n x 3] expressed in geographic
% east-north-up coordinates, into magnetic east-north-up coordinates.
%
% Seebany Datta-Barua
% 19 June 2013

function mENU = genu2menu(lat, lon, ht, yyyy, doy, gENU)

%% gENU to magnetic east, north, up (mENU).
[gmlat, gmlon] = geo2mag(lat, lon, yyyy, doy);
colat = (90-gmlat)*pi/180;
[sinD, cosD, Br, Btheta, Bphi] = trig_decl(colat, gmlon*pi/180, ht);

minusD = -asin(sinD);
R = rot_new(minusD*180/pi,3);
for i = 1:size(R,3)
    mENU(i,:) = R(:,:,i)*gENU(i,:)';%[ dmspvx(i); dmspvy(i); dmspvz(i)];
end

% mENU = mENU'; % Make the output vector [n x 3]