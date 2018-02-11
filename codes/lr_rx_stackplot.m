function [] = lr_rx_stackplot(lrdata, year, doy, rcvr_struct, op_path, plot_str)

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

if nargin == 0
    year = '2013';
    doy = '342';
    plot_str = 'sp';
    load([mat_path, 'lrdata_', year, '_', doy, '.mat']);
end

spdata = lrdata{1};
rcvr_op = lrdata{4};

tlim = lrdata{5};
splim = lrdata{6};
signal = lrdata{7};
spmask = lrdata{8};
s4data = lrdata{9};
s4mask = lrdata{10};
s4lim = lrdata{11};
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

dailymean = mean(lrdata(:, 2));
for rr = 1:size(rcvr_struct, 1)
    rcvr_name = rcvr_struct(rr, :);
    if ismember(rcvr_name, rcvr_op, 'rows')
        [~, rr_op] = ismember(rcvr_name, rcvr_op, 'rows');
        svdata = lrdata(lrdata(:, 4) == rr_op, 1:2);
        if strcmp(plot_str, 's4')
            %             find(svdata(:,2)>=datalim)
            svdata(svdata(:, 2) >= datalim, 2) = datalim;
        end
        if isempty(svdata)
            rectangle('Position', [tlim(1), 0.25 + rr - 1, diff(tlim), 0.5], 'FaceColor', [0.4, 0.4, 0.4]);
            text(tlim(1)+diff(tlim)/2, rr-1+0.5, 'Data not available', ...
                'HorizontalAlignment', 'center');
            hold on;
        end
        svdata = sortrows(svdata, 1);
        time = svdata(:, 1);
        data = svdata(:, 2);
        meandata(:, rr) = mean(data);
        stddata(:, rr) = std(data);
        megadata{:, rr} = data;
        %         [time_e,sigmaphi_e] = discont_proc(time,data,dt);
        %         [color] = rx_color(rcvr_name);
        colormap(jet);
        caxis(clim);
        scatter(time, data/datalim+rr-1, 15, data, 'filled');
        hold on;
        %mean after mask
        %         if strcmp(doy,'342')
        %             plot(tlim([1 end]),(0.6/datalim+rr-1)*ones(size(tlim([1 end]))),'color',[0.5 0.5 0.5],'LineWidth',1);
        %             text(tlim(end),0.6/pi+rr-1,'$th_{stat}$','HorizontalAlignment','right','VerticalAlignment','bottom');
        %             plot(tlim([1 end]),(dailymean/2/pi+rr-1)*ones(size(tlim([1 end]))),'color',[0.5 0.5 0.5],'LineWidth',1);
        %             text(tlim(end),dailymean/2/pi+rr-1,'$th_{dyn}$','HorizontalAlignment','right','VerticalAlignment','bottom');
        %             hold on;
        %         end
        
    elseif strcmp(rcvr_name, 'ASTRArx') == 1
        rectangle('Position', [tlim(1), 0.25 + rr - 1, diff(tlim), 0.5], 'FaceColor', [0.4, 0.4, 0.4]);
        text(tlim(1)+diff(tlim)/2, rr-1+0.5, 'Data not available', ...
            'HorizontalAlignment', 'center');
        hold on;
    else
        rectangle('Position', [tlim(1), 0.25 + rr - 1, diff(tlim), 0.5], 'FaceColor', [0.4, 0.4, 0.4]);
        text(tlim(1)+diff(tlim)/2, rr-1+0.5, 'Not operational', ...
            'HorizontalAlignment', 'center');
        hold on;
    end
    text(tlim(1), rr-0.5, {rcvr_name; char(sitenum_struct(rr, :))}, 'HorizontalAlignment', 'right');
end

% save('./Dropbox/statistics.mat','meandata','stddata','megadata');
% cb = colorbar('YTick',mask:0.5:datalim);
cb = colorbar;
set(cb, prop, clim);
set(get(cb, 'YLabel'), 'String', cblabel, 'interpreter', 'latex');
% if lrdata(end, 1) - lrdata(1, 1) <= 0.5
%     tlim = lrdata([1 end], 1);
% end
xlim(tlim);
datetick('x', 'HH', 'keeplimits');
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Ytick', 1:rr, 'YTickLabel', ' ');
xlabel('Universal Time [hr]');
switch plot_str
    case 'sp'
        title(['Low-rate ', title_str, ' stackplot for DOY:', num2str(doy, '%03i'), ...
            ' of year:', num2str(year)]);
    case 's4'
        title(['Low-rate ', title_str, ' stackplot for DOY:', num2str(doy, '%03i'), ...
            ' of year:', num2str(year)]);
end

op_path2 = ['/data1/public/webplots/cases_quicklook/', num2str(year), ...
    sep, num2str(doy, '%03i'), sep];
comm = ['mkdir -p ', op_path2];
system(comm);
try
    saveas(gca, [op_path2, signal, '_', plot_str, '_lr_rx_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
    saveas(gca, [op_path2, signal, '_', plot_str, '_lr_rx_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.eps'], 'epsc2');
    saveas(gca, [op_path, signal, '_', plot_str, '_lr_rx_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
catch
    saveas(gca, [op_path, signal, '_', plot_str, '_lr_rx_stackplot_', num2str(year), '_', ...
        num2str(doy, '%03i'), '_', num2str(mask), '.png'], 'png');
end
close;
toc
end
