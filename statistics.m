clear all;
close all;
load sd_2014.mat;
% load sd.mat;
% profile on
if verLessThan('matlab', '8.4.0')
    op_path = '/Users/yangsu/Dropbox/';
else
    op_path = '/Users/yangsu/Dropbox/';
end

%%
figure;
for jjj = 1:365
    RH_(jjj, :) = datenum(['2014-01-01 ', RH(jjj, :)], 'yyyy-mm-dd HHMM');
    SH_(jjj, :) = datenum(['2014-01-01 ', SH(jjj, :)], 'yyyy-mm-dd HHMM');
end
plot([1:365]+datenum([2014, 1, 1, 0, 0, 0])-1, RH_, [1:365]+datenum([2014, 1, 1, 0, 0, 0])-1, SH_);
datetick('y', 'HH:MM');
datetick('x', 'mmm');
xlabel('day of year');
legend('Rise time', 'Set time');
title('Solar terminator for year 2014');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
% saveas(gcf,[op_path,'solarterm.jpg'],'jpg');
close;
% return;

%%
figure;
truetime = MSP_days(:, 1) + MEGA_MSP(:, 1) - 1;
truetime = truetime - 9 / 24;
id = find(truetime >= datenum([2014, 1, 1, 0, 0, 0]));
truetime = truetime(id, :);

MEGA_MSP = MEGA_MSP(id, :);

ttv = datevec(truetime);
doy = datenum([ttv(:, 1:3), zeros(size(ttv, 1), 3)]) - datenum([2014, 1, 1, 0, 0, 0]);
toffset = truetime - doy;
% datevec(toffset([1 end],:))
se = datenum([2014, 3, 20, 0, 0, 0]);
ss = datenum([2014, 6, 21, 0, 0, 0]);
fe = datenum([2014, 9, 23, 0, 0, 0]);
ws = datenum([2014, 12, 21, 0, 0, 0]);
spring = find(truetime(:, 1) >= se & truetime(:, 1) < ss);
summer = find(truetime(:, 1) >= ss & truetime(:, 1) < fe);
fall = find(truetime(:, 1) >= fe & truetime(:, 1) < ws);
winter = find(truetime(:, 1) < se | truetime(:, 1) >= ws);
colormap(jet);
for iii = 1:4
    switch iii
        case 1
            sub = spring;
            str = {'Spring'; 'Mar 20 to June 20'};
        case 2
            sub = summer;
            str = {'Summer'; 'June 21 to Sep 22'};
        case 3
            sub = fall;
            str = {'Fall'; 'Sep 23 to Dec 21'};
        case 4
            sub = winter;
            str = {'Winter'; 'Before Mar 20 or after Dec 21'};
    end
    subplot(2, 2, iii);
    scatter(toffset(sub, 1), MEGA_MSP(sub, 2), [], MEGA_MSP(sub, 2), '.');
    axis([datenum([2014, 1, 1, 0, 0, 0]), datenum([2014, 1, 2, 0, 0, 0]), 0, 2 * pi]);
    datetick('x', 'HH', 'keeplimits');
    title(str);
    xlabel('Local Time of the day [HH]');
    grid on;
    box on;
    clim = [0, 2 * pi];
    caxis(clim);
    cb = colorbar;
    set(get(cb, 'YLabel'), 'String', '\sigma_\Phi value [Rad]');
    if verLessThan('matlab', '8.4.0')
        set(cb, 'Clim', clim);
    else
        set(cb, 'limits', clim);
    end
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
% saveas(gcf,[op_path,'hours.jpg'],'jpg');
close;
% return;

%%
if ~verLessThan('matlab', '8.4.0')
    %2015a
    for rcvr_num = 4:6
        figure;
        subplot(2, 1, 1);
        histogram(MSP_days(MSP_days(:, 2) == rcvr_num, 3), 'binmethod', 'auto', ...
            'normalization', 'count');
        ylabel('Number of \sigma_{\Phi}');
        xlabel('\sigma_{\Phi} [rad]');
        title(['\sigma_{\Phi} histogram for year 2014, operational receivers = ', num2str(rcvr_num)]);
        subplot(2, 1, 2);
        histogram(MSP_days(MSP_days(:, 2) == rcvr_num, 3), 'binmethod', 'auto', 'binedge', 0:0.1:1, ...
            'normalization', 'count');
        ylabel('Number of \sigma_{\Phi}');
        xlabel('\sigma_{\Phi} [rad]');
        title(['zoomed-in for year 2014, operational receivers = ', num2str(rcvr_num)]);
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
        %         saveas(gcf,[op_path,'histogram_rcvr',num2str(rcvr_num),'.eps'],'epsc2');
        close;
    end
end
