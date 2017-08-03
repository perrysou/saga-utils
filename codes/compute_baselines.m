function [avgenu,avglla,avgxyz] = compute_baselines(NAVDATA,NAVDATA_O,init_time,xtime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
UTCTIME = gps2utc(NAVDATA(:,1:2));
TIME = datenum(UTCTIME);
abstime = xtime/24/3600+init_time;

XYZECEF = NAVDATA(TIME>=abstime(1)&TIME<=abstime(end),3:5);
XYZECEF_O = NAVDATA_O(TIME>=abstime(1)&TIME<=abstime(end),3:5);

avgxyz = mean(XYZECEF,1);
avgxyz_o = mean(XYZECEF_O,1);

avglla = ecef2lla(avgxyz,1);
avglla = [avglla(:,1:2)*180/pi avglla(:,3)]; 

avglla_o = ecef2lla(avgxyz_o,1);
olat = avglla_o(:,1)*180/pi;
olon = avglla_o(:,2)*180/pi;
oalt = avglla_o(:,3);
avgenu = xyz2enu_new(avgxyz,olat,olon,oalt);
end

