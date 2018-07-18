function [] = lr_hr_sigmaphi_plot(tspan, prn, home_dir, cases_folder, rcvr_op, year, doy, sep, signal, signal_type)
% interval for each receiver
tspan_utc = datevec(tspan);
tstt_hr = tspan(1, :) - 50 / 24 / 3600;
tend_hr = tspan(2, :) + 50 / 24 / 3600;
tspan_utc_hr = datevec([tstt_hr; tend_hr]);
tstt_hr = tspan_utc_hr(1, :);
tend_hr = tspan_utc_hr(2, :);
tstt = tspan_utc(1, :);
tend = tspan_utc(2, :);
hour = tspan_utc(1, 4);
RCVRNAME_LR = [];
RCVRNAME_HR = [];
%pull out read low-rate scintdata over this time span for all receivers
for rr = 1:size(rcvr_op, 1)
    rcvr_name = rcvr_op(rr, :);
    %         if strcmp(cases_folder(end-4:end-1),'pfrr')
    %             %folder_path for 2013 Poker Flat data
    %             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
    %             in_path = strcat([cases_folder,rcvr_name,sep,year,sep,doy,sep]);
    %         else
    %             %folder_path for 2013 Calgary data
    %             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);
    %             in_path = strcat([cases_folder,rcvr_name,sep,doy,sep]);
    %         end
    [in_path, op_path] = inoutpath(cases_folder, home_dir, year, doy, rcvr_name);
    outfilename = strcat('prn_files_', signal, sep, 'fltrdata.mat');
    load([op_path, outfilename]);
    scintdata = SPM_M;
    [color] = rx_color(rcvr_name);
    if ~isempty(scintdata)
        svdata = scintdata(scintdata(:, 3) == prn & scintdata(:, 1) >= tspan(1) & scintdata(:, 1) <= tspan(2), :);
        if ~isempty(svdata)
            %plot low-rate sigmaphi
            subplot(3, 1, 1)
            init_time = datenum([tspan_utc(1, 1:3), hour, 0, 0]);
            plot((svdata(:, 1) - init_time)*24*3600, svdata(:, 2), 'Color', color);
            grid on;
            xlim((tspan - init_time)*24*3600);
            %         datetick('x',':MM:SS');
            title(['Low Rate standard deviation of phase for PRN', num2str(prn)]);
            xstring = strcat({'Time [s] after '}, num2str(hour), ...
                ':00:00 UT on: ', datestr(tspan_utc(1, :), 'mm/dd/yy'));
            xlabel(xstring);
            ylabel('Low Rate \sigma_\phi [rad]');
            hold on;
            RCVRNAME_LR = [RCVRNAME_LR; rcvr_name];
            legend(RCVRNAME_LR, 'Location', 'NorthEastOutside');
        end
    end
    
    
    % Obtain and save high rate filtered S4 & Sigmahi data for the 300s
    DATAM2 = Fn_ReadHighRate_CASESdata_sdb(prn, op_path, cases_folder, rcvr_name, signal_type, tstt_hr, tend_hr);
    set_plot = 'B';
    if ~isempty(DATAM2)
        Fn_Plot_HighRate_CASESdata_sdb(prn, op_path, signal_type, set_plot);
        load([op_path, 'hr_prn_files_L1CA/HR_Scintdata_PRN', num2str(prn), '.mat']);
        %plot high rate sigmaphi
        subplot(3, 1, 2)
        dt = 5;
        [time_e, data_e] = discont_proc(data_PRN(:, 1), data_PRN(:, 3), dt);
        maxsp(rr) = max(abs(data_e));
        plot(time_e-hour*3600, data_e, 'Color', color);
        xlim((tspan - init_time)*24*3600);
        str = strcat('High rate standard deviation of phase for signal: ', signal, ', PRN:', num2str(prn));
        title(str)
        xstring = strcat({'Time [s] after '}, num2str(hour), ...
            ':00:00 UT on: ', datestr(tspan_utc(1, :), 'mm/dd/yy'));
        xlabel(xstring);
        ylabel('High Rate \sigma_\phi [rad]');
        if max(maxsp) > 0.5 / 1.2
            ylim([0, max(maxsp) * 1.2]);
        else
            ylim([0, 0.5]);
        end
        grid on;
        hold on;
        RCVRNAME_HR = [RCVRNAME_HR; rcvr_name];
        legend(RCVRNAME_HR, 'Location', 'NorthEastOutside');
        
        %plot high rate S4
        subplot(3, 1, 3)
        dt = 5;
        [time_e, data_e] = discont_proc(data_PRN(:, 1), data_PRN(:, 2), dt);
        maxs4(rr) = max(abs(data_e));
        plot(time_e-hour*3600, data_e, 'Color', color);
        xlim((tspan - init_time)*24*3600);
        str = strcat('High rate S_4 of power for signal: ', signal, ', PRN:', num2str(prn));
        title(str)
        xstring = strcat({'Time [s] after '}, num2str(hour), ...
            ':00:00 UT on: ', datestr(tspan_utc(1, :), 'mm/dd/yy'));
        xlabel(xstring);
        ylabel('High Rate S_4');
        if max(maxs4) > 0.5 / 1.2
            ylim([0, max(maxs4) * 1.2]);
        else
            ylim([0, 0.5]);
        end
        grid on;
        hold on;
        legend(RCVRNAME_HR, 'Location', 'NorthEastOutside');
    end
end

%save the low-rate plot
str_name = strcat(signal, '_LR&HR_SPH_PRN', num2str(prn));
outdir = strcat(op_path);
plotfile = strcat(outdir, str_name, '.eps');
saveas(gcf, plotfile, 'epsc2');
close;
end