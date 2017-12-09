function [] = getTEC()
close all;
clear;
dbstop if error;
t0 = datenum([2017, 8, 21, 16, 30, 0]);
tf = datenum([2017, 8, 21, 20, 30, 0]);
label_tec = 'TECU [$10^{16}m^{-2}$]';
label_time = 'HH:MM [UTC]';
title_tec = '1Hz Equivalent Vertical TEC measurements';
textkw = {'verticalalignment', 'bottom'};
el_threshold = 30;
% t0 = datenum([2017, 8, 21, 18, 20, 0]);
% tf = datenum([2017, 8, 21, 18, 40, 0]);

%%
filename = 'dataout_170821_1642_1942.o';
% filename = 'dataout_2013_342.obs';
% filename = 'mojk2330.17o';
% filename = 'dataout_2013_342_0326.obs';
% filename = 'brew2330.17o';
% filename = 'dsrc2330.17o';

[obs_time_rinex, prn_rinex, tecl_rinex, tecp_rinex] = getRinex(filename);
tecl_rinex = [zeros(1, length(prn_rinex)); diff(tecl_rinex, 1, 1)];

%%
[~, ipath, ~] = ver_chk;
cd(ipath);
load('SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/ionodata.mat');
time_iono = datenum(gps2utc([IONODATA(:, 1), IONODATA(:, 2)]));
isValid_iono = time_iono >= t0 & time_iono <= tf; %& IONODATA(:, end) == prn;
t_iono = time_iono(isValid_iono);
IONODATA = IONODATA(isValid_iono, :);
[t_iono_unique, ifirst, ~] = unique(t_iono, 'stable');

atec_rinex0 = NaN * ones(length(t_iono_unique), length(prn_rinex));
rtec_rinex0 = NaN * ones(length(t_iono_unique), length(prn_rinex));
for iu = 1:length(t_iono_unique)
    if iu == length(t_iono_unique)
        ilast = length(t_iono);
    else
        ilast = ifirst(iu+1) - 1;
    end
    prn_seq = IONODATA(ifirst(iu):ilast, end);
    atec_seq = IONODATA(ifirst(iu):ilast, 3);
    rtec_seq = IONODATA(ifirst(iu):ilast, 4);
    for sv = 1:length(prn_seq)
        iprn = prn_rinex == prn_seq(sv);
        atec_rinex0(iu, iprn) = atec_seq(sv);
        rtec_rinex0(iu, iprn) = rtec_seq(sv);
    end
end
atec_rinex = NaN * ones(length(obs_time_rinex), length(prn_rinex));
rtec_rinex = NaN * ones(length(obs_time_rinex), length(prn_rinex));
for sv = 1:length(prn_rinex)
    t1 = t_iono_unique(~isnan(atec_rinex0(:, sv)));
    if ~isempty(t1)
        atec1 = atec_rinex0(~isnan(atec_rinex0(:, sv)), sv);
        rtec1 = rtec_rinex0(~isnan(rtec_rinex0(:, sv)), sv);
        atec_rinex(:, sv) = interp1(t1, atec1, obs_time_rinex);
        rtec_rinex(:, sv) = interp1(t1, rtec1, obs_time_rinex);
    end
end


load('SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/txinfodata.mat');
time_tx = datenum(gps2utc([DATA_el(:, 1), DATA_el(:, 2)]));
isValid_el = time_tx >= t0 & time_tx <= tf & DATA_el(:, 3) >= 0; % & DATA_el(:, 4) == prn
% el_tx = deg2rad(DATA_el(isValid_el, 3));
t_tx = time_tx(isValid_el);
[t_tx_unique, ifirst, ~] = unique(t_tx, 'stable');
el_rinex0 = NaN * ones(length(t_tx_unique), length(prn_rinex));
for iu = 1:length(t_tx_unique)
    if iu == length(t_tx_unique)
        ilast = length(t_tx);
    else
        ilast = ifirst(iu+1) - 1;
    end
    prn_seq = DATA_el(ifirst(iu):ilast, 4);
    el_seq = DATA_el(ifirst(iu):ilast, 3);
    for sv = 1:length(prn_seq)
        iprn = prn_rinex == prn_seq(sv);
        el_rinex0(iu, iprn) = el_seq(sv);
    end
end

el_rinex = NaN * ones(length(obs_time_rinex), length(prn_rinex));
for sv = 1:length(prn_rinex)
    t1 = t_tx_unique(~isnan(el_rinex0(:, sv)));
    el1 = el_rinex0(~isnan(el_rinex0(:, sv)), sv);
    el_rinex(:, sv) = interp1(t1, el1, obs_time_rinex);
end
%%
figure;
[tecveq_rinex, of_rinex] = stec2vtec(tecl_rinex, tecp_rinex, deg2rad(el_rinex), prn_rinex);
subplot(1, 1, 1);
plot(obs_time_rinex, of_rinex);
[vmax, imax] = max(of_rinex);
text(obs_time_rinex(imax), vmax, num2str(prn_rinex), textkw{:});
title('Oblique Factors as a Function of IPP Zenith $\chi$');
xlim([t0, tf]);
ylabel('$\cos \chi$');
xlabel(label_time);
datetick('x', 'HH:MM', 'keeplimits');
tightfig;
saveas(gcf, ['SolarEclipse/PRNall_of_1Hz.pdf']);

tecveq_rinex(el_rinex <= el_threshold) = NaN;
figure;
subplot(1, 1, 1);
plot(obs_time_rinex, tecveq_rinex);
[tecmax, imax] = max(tecveq_rinex);
text(obs_time_rinex(imax), tecmax, num2str(prn_rinex), textkw{:});
title(title_tec);
xlim([t0, tf]);
ylabel(label_tec);
xlabel(label_time);
datetick('x', 'HH:MM', 'keeplimits');
tightfig;
saveas(gcf, ['SolarEclipse/PRNall_rinex_TEC_1Hz.pdf']);

figure;
subplot(1, 1, 1);
[tecveq_, of_] = stec2vtec(rtec_rinex, atec_rinex, deg2rad(el_rinex), prn_rinex);
tecveq_(el_rinex <= el_threshold) = NaN;
plot(obs_time_rinex, tecveq_);
[tecmax, imax] = max(tecveq_);
text(obs_time_rinex(imax), tecmax, num2str(prn_rinex), textkw{:});
title(title_tec);
xlim([t0, tf]);
ylabel(label_tec);
xlabel(label_time);
datetick('x', 'HH:MM', 'keeplimits');
tightfig;
saveas(gcf, ['SolarEclipse/PRNall_saga_TEC_1Hz.pdf']);
% keyboard;
% close all;
%%
load('SolarEclipse/mat/grid108/2017/233/prn_files_L1CA/channeldata.mat');
isSlipped = CHANNELDATA(:, 6) == 1;
prns_slipped = unique(CHANNELDATA(isSlipped, end), 'stable');
% isSlipped = CHANNELDATA(:, 6) == NaN;

time_ch = datenum(gps2utc([CHANNELDATA(:, 1), CHANNELDATA(:, 2)]));

figure;
for iprn = 1:length(prn_rinex)
    isBetween = time_ch >= t0 & time_ch <= tf & CHANNELDATA(:, end) == prn_rinex(iprn);
    %     t_ch = time_ch(isBetween);
    rows_l1 = isBetween & CHANNELDATA(:, end-1) == 0; % CHANNELDATA(:, end) == prn
    rows_l2 = isBetween & CHANNELDATA(:, end-1) == 2;
    
    p1 = CHANNELDATA(rows_l1, 4); p2 = CHANNELDATA(rows_l2, 4);
    l1 = CHANNELDATA(rows_l1, 3); l2 = CHANNELDATA(rows_l2, 3); % cycles
    
    t_l1 = time_ch(rows_l1);
    t_l2 = time_ch(rows_l2);
    [t_common, i1, i2] = intersect(t_l1, t_l2);
    if isempty(t_common)
        continue;
    end
    [tecl, tecp] = getStec(l1(i1), l2(i2), p1(i1), p2(i2));
    tecl = [0; diff(tecl)];
    el_ = interp1(t_tx_unique, el_rinex0(:, iprn), t_common);
    [tecveq, of] = stec2vtec(tecl, tecp, deg2rad(el_), prn_rinex(iprn));
    
    %plot
    %     subplot(1, 1, 1);
    %     [t_l1_e, tecl_e] = discont_proc(t_common, tecl, 2 / 24 / 3600);
    %     plot(t_l1_e, [0; diff(tecl_e)]); hold on;
    %     title(['Unsmoothed 1Hz Carrier TEC measurements']);
    % %     xlim([t0, tf]);
    %     datetick('x', 'HH:MM', 'keeplimits');
    
    subplot(1, 1, 1);
    %     [t_l1_e, tecveq_e] = discont_proc(t_common, tecveq, 2 / 24 / 3600);
    tecveq(el_ <= el_threshold) = NaN;
    plot(t_common, tecveq); hold on;
    [tecmax, imax] = max(tecveq);
    text(t_common(imax), tecmax, num2str(prn_rinex(iprn)), textkw{:});
end
xlim([t0, tf]);
title(title_tec);
ylabel(label_tec);
xlabel(label_time);
datetick('x', 'HH:MM', 'keeplimits');
tightfig;
saveas(gcf, ['SolarEclipse/PRNall_raw_TEC_1Hz.pdf']);

return;
%%
% high-rate
load(['SolarEclipse/mat/grid108/2017/233/hr_prn_files_L1CA/PRN', num2str(prn), '_1800UT_zoom1.mat']);
time_l1 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l1_hr = DATAM2(time_l1 >= t0 & time_l1 <= tf, 4);

load(['SolarEclipse/mat/grid108/2017/233/hr_prn_files_L2CL/PRN', num2str(prn), '_1800UT_zoom1.mat']);
time_l2 = datenum(gps2utc([DATAM2(:, 2), DATAM2(:, 3)]));
l2_hr = DATAM2(time_l2 >= t0 & time_l2 <= tf, 4);

t_l1_hr = time_l1(time_l1 >= t0 & time_l1 <= tf);
t_l2_hr = time_l2(time_l2 >= t0 & time_l2 <= tf);

% figure;
plot(t_l1_hr, -l1_hr);
plot(t_l2_hr, -l2_hr);

[tecl_hr, ~] = getStec(-l1_hr, -l2_hr, NaN, NaN);
% [t_l1_hr_e, tecl_hr_e] = discont_proc(t_l1_hr, tecl_hr, 0.015 / 24 / 3600);

%hr
figure;
plot(t_l1_hr, tecl_hr, '.');
xlim([t0, tf]);
datetick('x', 'HH:MM', 'keeplimits');
title(['Unsmoothed 100Hz Carrier Phase TEC measurements, PRN', num2str(prn)]);
tightfig;
saveas(gcf, ['SolarEclipse/PRN', num2str(prn), '_TEC_100Hz.png']);
end