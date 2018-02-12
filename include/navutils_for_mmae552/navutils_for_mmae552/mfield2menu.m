% function v_menu = mfield2menu(lat, lon, yyyy, doy, mfield)
% rotates a vector mfield = [n x 3] expressed in magnetic field-aligned
% and perpendicular coordinates, into magnetic east-north-up coordinates.
% lat and lon are geographic coordinates.
%
% Seebany Datta-Barua
% 19 June 2013

function v_menu = mfield2menu(lat, lon, yyyy, doy, mfield)

%% mENU to field-perpendicular east, field-perpendicular north+up, 
% and field-parallel.
[gmlat, gmlon] = geo2mag(lat, lon, yyyy, doy);
colat = (90-gmlat)*pi/180;

[sinI, cosI] = trig_dip(colat);
% A positive rotation in this direction.
I = asin(sinI);
R = rot_new(I*180/pi,1);
for i = 1:size(R,3)
    v_menu(i,:) = R(:,:,i)*mfield(i,:)';
  
end
% v_mfield = v_mfield';