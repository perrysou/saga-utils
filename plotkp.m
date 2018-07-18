close all;
clear all;
op_path = '/Users/yangsu/Dropbox/research-misc/[MTSSP]/2015_Yang_MTSSP/figs/';
time = datenum([2015, 3, 17, 1.5, 0, 0]):1 / 24 * 3:datenum([2015, 3, 17, 22.5, 0, 0]);
%planetary K
kpdata = [2, 5, 6, 6, 8, 8, 7, 8];
%college K (high lat)
kcdata = [2, 3, 6, 6, 8, 8, 7, 5];
kpdata = kcdata;
weaki = find(kpdata < 4);
moderatei = find(kpdata == 4);
strongi = find(kpdata > 4);
hold on;
bar(time, kpdata, 1);
% b1 = bar(time(weaki),kpdata(weaki),'g');
% b2 = bar(time(moderatei),kpdata(moderatei),'y');
% b3 = bar(time(strongi),kpdata(strongi),'r');
% legend({'K<4','K=4','K>4'});
axis([time(1) - 1.5 / 24, time(end) + 1.5 / 24, 0, 9]);
grid on;
title('Estimated College K-index for High Latitude, updated every 3 hours')
xlabel('Universal Time [HH:MM], March 17, 2015');
ylabel('K-index value')
datetick('x', 'keeplimits');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
% saveas(gcf,[op_path,'Kp.eps'],'epsc2');
% close;
