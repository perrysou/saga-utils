function [tec] = getVerticalTEC(fArr, phiArr, rhoArr
t0 = datenum([2017 8 21 18 20 0]);
tf = datenum([2017 8 21 18 40 0]);

load('/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/channeldata.mat');

time_ch = datenum(gps2utc([CHANNELDATA(:, 1), CHANNELDATA(:, 2)]));
rows_l1 = find(CHANNELDATA(:, end)==12&CHANNELDATA(:, end-1)==0&time_ch>=t0&time_ch<=tf);
rows_l2 = find(CHANNELDATA(:, end)==12&CHANNELDATA(:, end-1)==2&time_ch>=t0&time_ch<=tf);

p1 = CHANNELDATA(rows_l1, 4);
l1 = CHANNELDATA(rows_l1, 3); % cycles
p2 = CHANNELDATA(rows_l2, 4);
l2 = CHANNELDATA(rows_l2, 3); % cycles
t_l1 = time_ch(rows_l1);
t_l2 = time_ch(rows_l2);
all(t_l1 == t_l2)

rhoArr = [p1 p2];
phiArr = [l1 l2];

[tecl, tecp] = getStec(l1, l2, p1, p2);

% high-rate
load('/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/hr_prn_files_L1CA/PRN12_1800UT_zoom1.mat');
time_l1 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l1_hr = DATAM2(time_l1>=t0&time_l1<=tf, 4);

load('/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/hr_prn_files_L2CL/PRN12_1800UT_zoom1.mat');
time_l2 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l2_hr = DATAM2(time_l2>=t0&time_l2<=tf, 4);

t_l1_hr = time_l1(time_l1>=t0&time_l1<=tf);
t_l2_hr = time_l2(time_l2>=t0&time_l2<=tf);
t_l1_hr == t_l2_hr;

[tecl_hr, ~] = getStec(l1_hr, l2_hr, NaN, NaN);

%hr
figure;
plot(t_l1_hr, tecl_hr);
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');
title('Unsmoothed 100Hz Carrier Phase TEC measurements, PRN12');
tightfig;
saveas(gcf, '/Users/yangsu/Dropbox/SolarEclipse/PRN12_TEC_100Hz.pdf');
%lr
figure;
subplot(2,1,1);
plot(t_l1, tecl);
title('Unsmoothed 1Hz Carrier Phase TEC measurements, PRN12');
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');

subplot(2,1,2);
plot(t_l1, tecp);
title('Unsmoothed 1Hz Code Pseudorange TEC measurements, , PRN12');
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');
tightfig;
saveas(gcf, '/Users/yangsu/Dropbox/SolarEclipse/PRN12_TEC_1Hz.pdf');

close all;
end