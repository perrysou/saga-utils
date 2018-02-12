% Script to check that v_ipp is small compared to ion drift velocity, for
% SAGA study.  I'll use the 12/8/13 case study, the txinfo file has az, el.
%
% From http://cases.astraspace.net/documentation/txinfodef.txt:
% ============================ txinfodef.txt ==================================
% This file defines the columns of data in the txinfo.log files produced by the
% GRID software receiver.  Each txinfo.log file contains time-stamped
% transmitter information in the form of azimuth angle, elevation angle, and
% health status for tracked transmitters. See channeldef.txt for a definition of
% ORT.  ORT time stamps indicate the time at which the transmitter information
% applies.
% =============================================================================
% 
% 
% Column        Quantity
% 
% 1 ----------- ORT week number. 
% 
% 2 ----------- ORT whole seconds of week.
% 
% 3 ----------- ORT fractional second. 
% 
% 4 ----------- Azimuth in degrees.
% 
% 5 ----------- Elevation in degrees.
% 
% 6 ----------- Health status, where zero is nominal.
% 
% 7 ----------- System:
%               0       GPS
%               1       GALILEO
%               2       GLONASS
%               3       CDMA
% 
% 8 ----------- Transmitter identification number (TXID).
% 
% =============================================================================
%
% Seebany Datta-Barua
% 9 Jun 2016

clear
r2d = 180/pi;
d2r = pi/180;

svinfo = load('/Users/seebany/Dropbox/Data/CASES/pfrr/grid108/2013/342/txt/txinfo_dataout_2013_342_0326.log');

prn = svinfo(:,8);

rows = find(prn == 23);

az = svinfo(rows,4)*pi/180; 
el = svinfo(rows,5)*pi/180; 

% Pick IIT-1 as ground station.  From Datta-Barua et al., 2015 GRL paper:
stnlat = 65.1265;%*pi/180;
stnlon = -147.4968;%*pi/180;
stnht = 0;

alt_H_iono = 350e3; %m
t = svinfo(rows,2)+svinfo(rows,3);
delT = mean(diff(t));

% Method 1: use Jiyun's function.
% Something is screwy here; the result seems to be 16 km/s!
% [ipplat ipplon V_ipp direction_ipp] = iono_as_getIPP(az, el, [stnlat, stnlon], H_iono, delT);

% Method 2: use my version of sill.m
[ipplat, ipplon] = sill(az, el, stnlat, stnlon, stnht, alt_H_iono); %output deg
xyz = wgslla2xyz(ipplat, ipplon, alt_H_iono*ones(numel(ipplat),1)); %output m

enu = xyz2enu(xyz,stnlat,stnlon,stnht);
denu = diff(enu);
dist = sqrt(denu(:,1).^2 + denu(:,2).^2);
v = dist./diff(t);

% Method 3: sort through what's going on in Jiyun's function.
for i = 1:numel(ipplat)-1
    enu2(i,:) = xyz2enu(xyz(i+1,:), ipplat(i), ipplon(i), 0);%alt_H_iono);
end
dist2 = sqrt(enu2(:,1).^2 + enu2(:,2).^2);
vbad = dist2/delT;

% So Methods 2 and 3 are now the same.

IPPlla2 = [ipplat(2:end)*pi/180 ipplon(2:end)*pi/180 alt_H_iono*ones(numel(ipplat)-1,1)];
IPPlla1 = [ipplat(1:end-1)*pi/180 ipplon(1:end-1)*pi/180 alt_H_iono*ones(numel(ipplat)-1,1)];
% 
% Next is same as xyz.
IPPecef2 = lla2ecef(IPPlla2);

% For some reason, using IPPlla1 in the next line gives a weird north
% motion.  Also the down direction is too large a magnitude.
% I Suspect the ecef2ned function I'm using is not the one it expects or
% something.  This would be the same problem to happen if using Jiyun's
% function on my same computer, with the same ecef2ned I'm using now.
    IPPned = ecef2ned(IPPecef2, IPPlla1);%[stnlat, stnlon, stnht]*pi/180);%
% Next line is not correct if using a fixed stnlat/lon/ht; need to compute
% a change in position.
%     
%     dist = sqrt((IPPned(:,1)).^2+(IPPned(:,2)).^2);
%     v = dist/delT;
%     direction=acos(IPPned(:,2)./dist)*r2d;
