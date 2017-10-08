function [AZb, ELb, beamid, op_path] = plotPFISR_NeTe(t1, t2, flag)
t0 = [2013, 12, 8, 3, 30, 0];
tf = [2013, 12, 8, 4, 30, 0];
flag = 'Ne';
[~, op_path, mat_path] = ver_chk;
try
    load([mat_path, 'PFISR_data2.mat'], 'PFISR_data');
catch
    PFISR_data = importdata('PFISR_data2.txt');
    save([mat_path, 'PFISR_data2.mat'], 'PFISR_data');
end

% columns
% 1:DATENUM,2:AZM,3:ELM,4:BEAMID,5:RANGE,6:NEL,7:DNEL,8:TI,9:DTI,10:TE,11:DTE

% set up flags
if strcmp(flag, 'Ne')
    cols = 6;
end

delta = 0;
if delta
    deltaflag = 'delta';
else
    deltaflag = 'normal';
end

%creat folders for figures
op_path = strjoin({op_path, 'research-misc', 'ION2017', 'figs', ...
    flag, deltaflag, filesep}, filesep);
mkdir = strjoin({'mkdir -p', op_path});
system(mkdir);

data = [datenum(PFISR_data(:, 1:6)), PFISR_data(:, 7:end)];
data = data(data(:, 1) <= datenum(tf)+360/24/3600 & data(:, 1) >= datenum(t0)-360/24/3600,:);
if ~isempty(data)
    beamid = unique(data(:, 4), 'stable');
    for ibeam = length(beamid)
        AZb(ibeam,:) = unique(data(data(:, 4) == beamid(ibeam), 2));
        ELb(ibeam,:) = unique(data(data(:, 4) == beamid(ibeam), 3));
        beamstr = strjoin({num2str(beamid(ibeam)), ...
            [num2str(AZb(ibeam,:)), '$^\circ$ az'], [num2str(ELb(ibeam,:)), '$^\circ$ el']}, ', ');
%         subplot(2,1,2);
        data_beam = data(data(:, 4) == beamid(ibeam),:);
        ranges = unique(data_beam(:, 5), 'stable');
        times = unique(data_beam(:, 1), 'stable');
        ne = 10.^data_beam(:, cols);
        te = data_beam(:, 10);
        negrid = reshape(ne, length(ranges), []);
%         continue;
        difftimenegrid = [diff(negrid, 1, 2), NaN(length(ranges), 1)] ./ ...
            [diff(times) * 24 * 3600; NaN]';
        diffrangenegrid = [diff(negrid, 1, 1); NaN(1, length(times))] ./ ...
            [diff(ranges) * 10e3; NaN];
        tegrid = reshape(te, length(ranges), []);
        [timegrid, rangegrid] = meshgrid(times, ranges);
        if delta
            pcolor(timegrid, rangegrid, difftimenegrid);
            caxis([-1.5 * 10 ^ 9, 1.5 * 10 ^ 9]);
            title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                'PFISR Electron Density Variation $\Delta {N_e} / \Delta t, \Delta t \approx 180 s$}']);
            cb = colorbar;
            set(get(cb, 'YLabel'), 'String', ...
                '$\frac{\Delta N_e}{\Delta t} [m^{-3} s^{-1}]$', ...
                'interpreter', 'latex');
        else
            pcolor(timegrid, rangegrid, log10(negrid));
            caxis([10, 12]);
            title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                'PFISR Electron Density $N_e$}']);
            cb = colorbar;
            set(get(cb, 'YLabel'), 'String', ...
                '$\log_{10}{N_e} [m^{-3}]$', ...
                'interpreter', 'latex');
        end
        shading flat;       
        set(gca, 'layer', 'top');
        grid off;
        xlabel('Time [HH:MM UT]');
        ylabel('Altitude [km]');
        ylim([100, 600]);
        xlim(datenum([t1; t2]));
        datetick('x', 'HH:MM', 'keeplimits');
        plotname = strjoin({num2str(beamid(ibeam)), flag, deltaflag}, '_');
        plotpath = [op_path, plotname, '.png'];
        saveas(gcf, plotpath, 'png');
    end
else
    beamid = NaN;
    AZb = NaN;
    ELb = NaN;
end
% close all;
end


