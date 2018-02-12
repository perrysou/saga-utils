% function v_mfield = menu2mfield(lat, lon, yyyy, doy, mENU)
% rotates a vector gENU = [n x 3] expressed in geographic
% east-north-up coordinates, into magnetic east-north-up coordinates.
%
% Seebany Datta-Barua
% 19 June 2013

function v_mfield = menu2mfield(lat, lon, yyyy, doy, mENU)

%% mENU to field-perpendicular east, field-perpendicular north+up, 
% and field-parallel.
[gmlat, gmlon] = geo2mag(lat, lon, yyyy, doy);
colat = (90-gmlat)*pi/180;

[sinI, cosI] = trig_dip(colat);
minusI = -asin(sinI);
R = rot_new(minusI*180/pi,1);
for i = 1:size(R,3)
    v_mfield(i,:) = R(:,:,i)*mENU(i,:)';
  
end
% v_mfield = v_mfield';