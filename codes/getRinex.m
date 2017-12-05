dbstop if error;
clear;
[~, ipath, opath] = ver_chk;
cd(strjoin({ipath, 'SolarEclipse', 'CORS_RINEX_MOJK'}, filesep));
filename = 'dataout_170821_1642_1942.o';
% filename = 'mojk2330.17o';
% filename = 'dataout_2013_342_0326.obs';
% filename = 'brew2330.17o';
% filename = 'dsrc2330.17o';
outfmt = '.mat';
fpath = strjoin({ipath, 'SolarEclipse', 'CORS_RINEX_MOJK', filename}, filesep);
matpath = strjoin({ipath, 'SolarEclipse', 'CORS_RINEX_MOJK', [filename, outfmt]}, filesep);
if ~exist(matpath, 'file')
    [obs,obs_qual,type_str,clk_off] = rd_rnx_o6(filename);
else
    load(matpath);
end
[gps_time, prn, obsout] = refrinex3(obs, obs_qual, clk_off);
obs_time = datenum(gps2utc(gps_time(:,1:2), gps_time(:,3)));


for i = 1:size(type_str, 1)
    type = type_str(i, :);
    switch type
        case 'L1'
            l1 = obsout{i};
        case 'L2'
            l2 = obsout{i}; 
        case 'C1'
            p1 = obsout{i}; 
        case 'C2'
            p2 = obsout{i};
    end
end

ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
hold(ax1, 'on');
hold(ax2, 'on');
for j = 1:length(prn)       
    [tecl_, tecp_] = getStec(l1(:, j), l2(:, j), p1(:, j), p2(:, j));    
    plot(ax1, obs_time, tecl_, '.'); 
    plot(ax2, obs_time, tecp_, '.'); 
end
title(ax1, ['Phase Slant TEC, PRN:']);
title(ax2, ['Code Slant TEC, PRN:']);
legend(ax2, num2str(prn));
datetick(ax1, 'x', 'HH:MM');
datetick(ax2, 'x', 'HH:MM');
    