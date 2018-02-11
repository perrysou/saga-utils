function [] = lr_prn_stackplot(lrdata, year, doy, rcvr_struct, op_path, plot_str)

[prop, op_path, mat_path] = ver_chk;

% %where cases data are stored
% case_folder = '/data1/public/cases/pfrr/';
% %prefix of cases receiver ID
% rcvr_prefix = 'grid';
% %year
% year = '2013';
% %day of year
% doy = '342';
% %where output data and plots are stored
% home_dir = '/data1/home/ysu27/';
%
% if strcmp(case_folder(end-4:end-1),'pfrr')
%     %folder_path for 2013 Poker Flat data
%     op_path = [home_dir,'PFRR_Data/'];
% else
%     %folder_path for 2013 Calgary data
%     op_path = [home_dir,'Calgary_Data/'];
% end
%
% %load lrdata
% load([op_path,'lrdata_',year,'_',doy,'.mat']);
tic
sep = filesep;

% load([mat_path,'lrdata_',year,'_',doy,'.mat']);

spdata = lrdata{1};
rcvr_op = lrdata{4};
tlim = lrdata{5};
splim = lrdata{6};
signal = lrdata{7};
spmask = lrdata{8};
s4data = lrdata{9};
s4mask = lrdata{10};
s4lim = lrdata{11};
prnlist = (1:32)';
dt = 1 / 24 / 3600 * 100;

switch plot_str
    case 'sp'
        lrdata = spdata;
        mask = spmask;
        clim = [0, 2 * pi];
        datalim = clim(2);
        cblabel = '$\sigma_{\Phi}$ [rad]';
        mean_str = '$\sigma_{\Phi_{mean}} = $';
        title_str = '$\sigma_{\Phi}$';
    case 's4'
        lrdata = s4data;
        mask = s4mask;
        clim = [0, 0.5];
        datalim = clim(2);
        cblabel = '$S_4$';
        mean_str = '$S_{4_{mean}} = $';
        title_str = '$S_4$';
end

tlim1 = datenum([str2double(year), 0, 0, 0, 0, 0]) + str2double(doy);
tlim = [tlim1, tlim1 + 1];
sitenum_struct = rx2site(rcvr_struct);

for kk = 1:size(prnlist, 1)
    prn = prnlist(kk, :);
    svdata = lrdata(lrdata(:, 3) == prn, 1:2);
    if ~isempty(svdata)
        time = svdata(:, 1);
        data = svdata(:, 2);
        %     [time_e,sigmaphi_e] = discont_proc(time,data,dt);
        %             f1 = plot(time_e,sigmaphi_e/2/pi+kk-0.5,'LineWidth',1,'Color',color);
        hold on;
        colormap(jet);
        caxis(clim);
        scatter(time, data/datalim+kk-0.5, 13, data, 'filled');
        tmin = min(svdata(:, 1));
    else
        tmin = tlim(1);
    end
    text(tmin, kk, num2str(prn), 'VerticalAlignment', 'Baseline', 'HorizontalAlignment', ...
        'Right');
end
cb = colorbar;
set(cb, prop, clim);
set(get(cb, 'YLabel'), 'String', cblabel, 'interpreter', 'latex');
% if lrdata(end, 1) - lrdata(1, 1) <= 0.5
%     tlim = lrdata([1 end], 1);
% end
xlim(tlim);
datetick('x', 'HH', 'keeplimits');
box on;
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Ytick', 0.5:kk-0.5, 'YTickLabel', ' ');
xlabel('Universal Time [hr]');
switch plot_str
    case 'sp'
        title(['Low-rate ', title_str, ' prn-stackplot for DOY:', num2str(doy, '%03i'), ...
            ' of year:', num2str(year)]);
    case 's4'
        title(['Low-rate ', title_str, ' prn-stackplot for DOY:', num2str(doy, '%03i'), ...
            ' of year:', num2str(year)]);
end

op_path2 = ['/data1/public/webplots/cases_quicklook/', num2str(year), ...
    sep, num2str(doy, '%03i'), sep];
comm = ['mkdir -p ', op_path2];
system(comm);
try
    saveas(gca, [op_path2, signal, '_', plot_str, '_lr_prn_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
    saveas(gca, [op_path2, signal, '_', plot_str, '_lr_prn_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.eps'], 'epsc2');
    saveas(gca, [op_path, signal, '_', plot_str, '_lr_prn_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
catch
    saveas(gca, [op_path, signal, '_', plot_str, '_lr_prn_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
end
close;
toc
end
