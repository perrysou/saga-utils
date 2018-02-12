% Function [lat_sm, lon_sm] = ll2sm(abs_time,lat_geo,lon_geo)
% geodetic coordinates of a point to the solar magnetic frame.
% abs_time is absolute GPS time in seconds.  lat_geo and lon_geo in degrees.
%
% Seebany Datta-Barua
% 10.07.01

function [lat_sm, lon_sm] = ll2sm(abs_time, lat_geo, lon_geo)

D2R = pi/180;

% Greenwich hour angle.  Need to figure out what fmod and M_PI are.
gha = (mod(abs_time/3600.0,24.0) - 12.0) * pi/12.0

lat_se = lat_geo*D2R;
lon_se = lon_geo*D2R + gha;

% Location of magnetic pole taken from /vobs/tms/src/coord_xform.c at SE2SM
% and from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MAG_POLE_LAT = 78.6*D2R;
MAG_POLE_LON = -70.1*D2R;
MAG_INCLINATION = 11.5*D2R;

cos_pt_SE_lat = cos(lat_se);
cos_pt_SE_lon = cos(lon_se);
sin_pt_SE_lat = sin(lat_se);
sin_pt_SE_lon = sin(lon_se);

cos_pole_SE_lat = cos(MAG_POLE_LAT);
cos_pole_SE_lon = cos(MAG_POLE_LON + gha);
sin_pole_SE_lon = sin(MAG_POLE_LON + gha);

e_hat=[-cos_pole_SE_lat*sin_pole_SE_lon; cos_pole_SE_lat*cos_pole_SE_lon; 0.0];

e_hat_length = sqrt(e_hat(1)^2 + e_hat(2)^2);
e_hat = e_hat/e_hat_length;

e_hatx = [0 -e_hat(3) e_hat(2); e_hat(3) 0 -e_hat(1); -e_hat(2) e_hat(1) 0];

e_e1 = e_hat(1)*e_hat;

e_e2 = e_hat(2)*e_hat;

e_e3 = e_hat(3)*e_hat;

e_e = [e_e1 e_e2 e_e3];

T_SM_SE = eye(3) - sin(MAG_INCLINATION)*e_hatx + (1-cos(MAG_INCLINATION))*(e_e-eye(3));

r_SE = [cos_pt_SE_lat*cos_pt_SE_lon; cos_pt_SE_lat*sin_pt_SE_lon; sin_pt_SE_lat];

r_SM = T_SM_SE * r_SE;

length_r_SM = sqrt(r_SM(1)^2+r_SM(2)^2 + r_SM(3)^2);

lon_sm = atan2(r_SM(2),r_SM(1))/D2R;
lat_sm = asin(r_SM(3)/length_r_SM)/D2R;

return
