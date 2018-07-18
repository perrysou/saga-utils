function [AZ, EL, beamid] = readHDF5_Ne_Te_beam11()
%
% clear all;
% close all;

op_path = '/home/yang/Dropbox/';
%read data in datavec
tstt = [2013, 12, 8, 2, 0, 0];
tend = [2013, 12, 8, 6, 0, 0];
plot_str_v = ['Ne'; 'Te'];
%HDF5 filename
hdf5file = 'pfa131207.004.hdf5';
% h5disp(hdf5file);
info = h5info(hdf5file);
beam_struct = info.Groups(1).Groups.Groups;
for i = 11:size(beam_struct, 1)
    beamstr(i, :) = beam_struct(i).Name;
    beamid(i, :) = beamstr(i, end-5:end-1);
    beam_struct(i).Datasets(3).Name;
    tsec = h5read(hdf5file, [beamstr(i, :), '/', beam_struct(i).Datasets(3).Name]);
    tdt = datenum([1970, 1, 1, 0, 0, 0]) + double(tsec) / 24 / 3600;
    %     tvec = datevec([tdt(1);tdt(end)])
    ind = find(tdt <= datenum(tend)+10/24/3600 & tdt >= datenum(tstt)-10/24/3600);
    %     ind = [1:length(tdt)]';
    rc_1D = ind(end, :) - ind(1, :) + 1;
    topleft = [ind(1, :), 1];
    rc_2D = [ind(end, :) - ind(1, :) + 1, inf];
    str_1D = beam_struct(i).Groups(1).Name;
    str_2D = beam_struct(i).Groups(2).Name;
    az = h5read(hdf5file, [str_1D, '/azm'], ind(1, :), rc_1D);
    AZ(i, :) = unique(az, 'rows');
    el = h5read(hdf5file, [str_1D, '/elm'], ind(1, :), rc_1D);
    EL(i, :) = unique(el, 'rows');
    nel = h5read(hdf5file, [str_2D, '/nel'], topleft, rc_2D);
    te = h5read(hdf5file, [str_2D, '/te'], topleft, rc_2D);
    gdalt = h5read(hdf5file, [str_2D, '/gdalt'], topleft, rc_2D);
    gdalt = unique(gdalt, 'rows');
    
    for j = 1:size(plot_str_v, 1)
        plot_str = plot_str_v(j, :)
        %plot Ne/Te
        subplot(size(plot_str_v, 1), 1, j);
        set(gca, 'Layer', 'top');
        colormap(jet);
        switch plot_str
            case 'Ne'
                data = nel;
                clim = [10, 12];
                cbtick = 10:0.4:12;
                cbticklabel = num2str(cbtick', '%.1f');
                label_str = 'log_{10}Ne (m^{-3})';
            case 'Te'
                data = te;
                clim = [0, 3000];
                cbtick = 0:600:3000;
                cbticklabel = cbtick;
                label_str = 'Te (K)';
        end
        caxis(clim);
        for ii = 1:length(ind)
            for jj = 1:length(gdalt)
                if ii == length(ind)
                    xdata = [tdt(ind(ii)); tdt(ind(ii)) + 180 / 24 / 3600; tdt(ind(ii)) + 180 / 24 / 3600; tdt(ind(ii))];
                else
                    xdata = [tdt(ind(ii)); tdt(ind(ii+1)); tdt(ind(ii+1)); tdt(ind(ii))];
                end
                if jj == length(gdalt)
                    ydata = [gdalt(jj); gdalt(jj); gdalt(jj) + 22; gdalt(jj) + 22];
                else
                    ydata = [gdalt(jj); gdalt(jj); gdalt(jj+1); gdalt(jj+1)];
                end
                zdata = ones(4, 1);
                
                patch(xdata, ydata, zdata, 'EdgeColor', 'none', 'FaceColor', 'flat', 'CData', data(ii, jj), 'CDataMapping', 'scaled');
                hold on;
            end
        end
        title(['Beam:', num2str(i), ' (', num2str(AZ(i, :), '%04.1f'), '\circ az, ', num2str(EL(i, :), '%04.1f'), '\circ el)']);
        %     grid on;
        box on;
        set(gca, 'YTick', 150:50:500);
        ylabel('Altitude [km]');
        tsn = datenum(tstt);
        ten = datenum(tend);
        xticks = tsn:0.5 / 24:ten;
        set(gca, 'XTick', xticks);
        %     datevec(get(gca,'Xtick'))
        datetick('x', 'HH:MM', 'keepticks');
        axis([tsn, ten + 180 / 24 / 3600, 110, 560]);
        %     datevec(get(gca,'Xtick'))
        %         set(gca,'FontSize',14);
        if i == 11
            cb = colorbar('YTick', cbtick, 'YTickLabel', cbticklabel, 'CLim', clim);
            set(get(cb, 'YLabel'), 'String', label_str);
            xlabel('Time [HH:MM UT], DOY:342, YEAR:2013');
        end
    end
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
plotname = ['2013_342_Ne_Te_', ...
    datestr(tdt(ind(1)), 'HHMM'), '-', datestr(tdt(ind(end)), 'HHMM_UT')];
plotpath = [op_path, plotname, '.eps'];
saveas(gcf, plotpath, 'epsc2');
close;

% skyPlot([AZ AZ],[EL EL],str2double(cellstr(beamid)),'auto');
end
