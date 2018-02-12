% function gENU = menu2genu(lat, lon, ht, yyyy, doy, mENU)
% rotates a vector mENU = [n x 3] expressed in magnetic
% east-north-up coordinates, into geographic east-north-up coordinates.
%
% Seebany Datta-Barua
% 19 June 2013

function gENU = menu2genu(lat, lon, ht, yyyy, doy, mENU)

%% gENU to magnetic east, north, up (mENU).
[gmlat, gmlon] = geo2mag(lat, lon, yyyy, doy);
colat = (90-gmlat)*pi/180;
[sinD, cosD, Br, Btheta, Bphi] = trig_decl(colat, gmlon*pi/180, ht);

D = asin(sinD);
R = rot_new(D*180/pi,3);
for i = 1:size(R,3)
    gENU(i,:) = R(:,:,i)*mENU(i,:)';%[ dmspvx(i); dmspvy(i); dmspvz(i)];
end

% mENU = mENU'; % Make the output vector [n x 3]