function [] = getVerticalTEC()
t0 = datenum([2017, 8, 21, 18, 20, 0]);
tf = datenum([2017, 8, 21, 18, 40, 0]);

prn = 12;

load('/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/txinfodata.mat');
time_tx = datenum(gps2utc([DATA_el(:, 1), DATA_el(:, 2)]));
isValid_el = time_tx >= t0 & time_tx <= tf & DATA_el(:, 3) >= 30;
el_tx = DATA_el(isValid_el, 3);
t_tx = time_tx(isValid_el);
plot(t_tx, el_tx);

RE = 6.3781 * 10^6;
h = 350 * 1000;
chi = asin(RE / (RE + h) * sin(pi / 2 + deg2rad(el_tx)));
slant2vertical = cos(chi);


load('/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/channeldata.mat');

isSlipped = CHANNELDATA(:, 6) == 1;

time_ch = datenum(gps2utc([CHANNELDATA(:, 1), CHANNELDATA(:, 2)]));
isBetween = time_ch >= t0 & time_ch <= tf;

rows_l1 = isBetween & ~isSlipped & CHANNELDATA(:, end) == prn & CHANNELDATA(:, end-1) == 0;
rows_l2 = isBetween & ~isSlipped & CHANNELDATA(:, end) == prn & CHANNELDATA(:, end-1) == 2;

p1 = CHANNELDATA(rows_l1, 4);
l1 = CHANNELDATA(rows_l1, 3); % cycles
p2 = CHANNELDATA(rows_l2, 4);
l2 = CHANNELDATA(rows_l2, 3); % cycles
t_l1 = time_ch(rows_l1);
t_l2 = time_ch(rows_l2);

% figure;
% plot(t_l1); hold on; plot(t_l2);

% all(t_l1 == t_l2)
[tecl, tecp] = getStec(l1, l2, p1, p2);

tecl_smoothed = tecl - mean(tecl) + mean(tecp);

% high-rate
load(['/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/hr_prn_files_L1CA/PRN',num2str(prn), '_1800UT_zoom1.mat']);
time_l1 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l1_hr = DATAM2(time_l1 >= t0 & time_l1 <= tf, 4) / 2 / pi;

load(['/Users/yangsu/Dropbox/SolarEclipse/mat/grid108/2017/233/hr_prn_files_L2CL/PRN',num2str(prn), '_1800UT_zoom1.mat']);
time_l2 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l2_hr = DATAM2(time_l2 >= t0 & time_l2 <= tf, 4) / 2 / pi;

t_l1_hr = time_l1(time_l1 >= t0 & time_l1 <= tf);
t_l2_hr = time_l2(time_l2 >= t0 & time_l2 <= tf);

% figure;
% plot(t_l1_hr); hold on; plot(t_l2_hr);

[m, im] = nearest(t_l1_hr(1:10), t_l2_hr(1:11));
if length(t_l1_hr) < length(t_l2_hr)
    t_l2_hr_new = t_l2_hr((length(t_l2_hr) - length(t_l1_hr) + 1):end);
    l2_hr_new = l2_hr((length(t_l2_hr) - length(t_l1_hr) + 1):end);
    l2_hr = l2_hr_new;
end

[tecl_hr, ~] = getStec(l1_hr, l2_hr, NaN, NaN);
tecl_hr = tecl_hr - mean(tecl_hr) + mean(interp1(t_l1, tecp, t_l1_hr, 'linear', 'extrap'));
% [t_l1_hr_e, tecl_hr_e] = discont_proc(t_l1_hr, tecl_hr, 0.015 / 24 / 3600);

%hr
figure;
plot(t_l1_hr, tecl_hr, '.');
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');
title(['Unsmoothed 100Hz Carrier Phase TEC measurements, PRN', num2str(prn)]);
tightfig;
saveas(gcf, ['/Users/yangsu/Dropbox/SolarEclipse/PRN', num2str(prn), '_TEC_100Hz.pdf']);
%lr
figure;
subplot(3, 1, 1);
[t_l1_e, tecl_e] = discont_proc(t_l1, tecl, 2 / 24 / 3600);
plot(t_l1_e, tecl_e);
title(['Unsmoothed 1Hz Carrier Phase TEC measurements, PRN', num2str(prn)]);
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');

subplot(3, 1, 2);
[t_l1_e, tecp_e] = discont_proc(t_l1, tecp, 2 / 24 / 3600);
plot(t_l1_e, tecp_e);
title(['Unsmoothed 1Hz Code Pseudorange TEC measurements, PRN', num2str(prn)]);
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');

subplot(3, 1, 3);
[t_l1_e, tecl_smoothed_e] = discont_proc(t_l1, tecl_smoothed, 2 / 24 / 3600);
plot(t_l1_e, tecl_smoothed);
title(['Smoothed 1Hz Carrier Phase TEC measurements, PRN', num2str(prn)]);
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');

tightfig;
saveas(gcf, ['/Users/yangsu/Dropbox/SolarEclipse/PRN', num2str(prn), '_TEC_1Hz.pdf']);

close all;
end