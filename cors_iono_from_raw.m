% function [waas_time, Iphi_unbiased] = grace_iono_from_raw(utc0, minutespan, stn)
% computes
% the code and carrier measurements of ionosphere, levels the code to the
% carrier, removes the broadcast satellite biases, and removes the receiver
% bias by setting the minimum measurement to 0.  Adapted from
% grace_iono_from_raw.m.
%
% Seebany Datta-Barua
% 28 Mar 2007
% 8 Aug 2013 Adapting for use on raw WAAS data on IIT MacOSX. SDB.

function [waas_time, Ipr, Iccd, Iphi_unbiased] = cors_iono_from_raw(utc0, minutespan, stn)

load gpsconst lambda_1 lambda_2 gamma
yy = num2str(utc0(1));
mm = num2str(utc0(2));
dd = num2str(utc0(3));
if utc0(3) < 10
    dd = ['0', dd];
end
hr = num2str(utc0(4));
if utc0(4) < 10
    hr = ['0', hr];
end
minute = num2str(utc0(5));
sec = num2str(utc0(6));
clear utc0
utc0.year = str2num(yy);
utc0.mon = str2num(mm);
utc0.day = str2num(dd);
utc0.hour = str2num(hr);
utc0.min = str2num(minute);
utc0.sec = str2num(sec);
[GPS, ABS, RINEX] = convert_utc(utc0);
t0 = sec2waas(ABS.sec);

% Load the L1,L2,P1,P2 data.
% load(['C:\Documents and Settings\seebany\My Documents\data\CORS\'...
%     stn '_' yy(3:4) mm dd '.mat']);
load(['/Users/seebany/Data/supertruth/raw/121009/', ...
    stn, '_', yy(3:4), mm, dd, '.mat']);

waas_time = gps2waas(gps_time(:, 1:2));
% Limit data to timespan specified above.
timerows = find(waas_time >= t0-minutespan*60 & ...
    waas_time <= t0+minutespan*60);
waas_time = waas_time(timerows);
gps_time = gps_time(timerows, 1:2);

% Compute iono delays.
Ipr = (P2(timerows, :) - C1(timerows, :)) / (gamma - 1);
Iphi = (L1(timerows, :) * lambda_1 - L2(timerows, :) * lambda_2) / (gamma - 1);
Iccd = (C1(timerows, :) - L1(timerows, :) * lambda_1) / 2;

% May not need this part anymore.
prnlist = [10, 13, 23]; %prn;%[2 3 8 11 13 27 28 31];%[1 2 3 7 8 11 13 27 28 31];%[2 3 8 11 13 27 28 31];%
[c, iprn, ib] = intersect(prn, prnlist);
for i = 1:length(iprn)
    Iphi(:, iprn(i)) = detectslips(gps_time, Iphi(:, iprn(i)));
    
    % Level each continuous segment to remove integer ambiguity.
    [startrows, endrows] = find_seg(Iphi(:, iprn(i)));
    
    for j = 1:length(startrows)
        %         [startrows(j) endrows(j)]
        Iphi(startrows(j):endrows(j), iprn(i)) = ...
            Iphi(startrows(j):endrows(j), iprn(i)) - ...
            mean(Iphi(startrows(j):endrows(j), iprn(i)));
        Iphi(startrows(j):endrows(j), iprn(i)) = ...
            Iphi(startrows(j):endrows(j), iprn(i)) + ...
            mean(Ipr(startrows(j):endrows(j), iprn(i)));
        % plot((waas_time-t0)/60, Iphi(:,iprn(i)),'g.-')
    end
end
% % Apply broadcast Tgd from Rinex data from that day.
% Iphi_no_tgd = rm_tgd(Iphi, prn, num2str(RINEX.day), yy(end-1:end));
% Ipr_no_tgd = rm_tgd(Ipr, prn, num2str(RINEX.day), yy(end-1:end));

% % Estimate rx IFB only with the reliable prns.
% bias = est_rx_bias(stn);
% Iphi_unbiased = Iphi_no_tgd - bias;
% Ipr_unbiased = Ipr_no_tgd - bias;
Iphi_unbiased = Iphi;