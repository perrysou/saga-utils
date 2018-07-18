function MEGA_LZ = spectral(xdata, tpeak, rcvr_op, sitenum_op, combos, ...
    AZ, ZE, vmag, vang, doy, year, tslist, telist, tt, prn, init_time)
warning off;
dbstop if error;
clear all;
load lz.mat;

[~, op_path, ~] = ver_chk;
% v_ccmin = [0.6];
% v_dtau = 60;

% flags used in the analysis
normalize = 0;
if normalize
    normflag = 'normalized';
else
    normflag = 'non-normalized';
end
fitall = 0;
if fitall
    fitflag = 'fitall';
else
    fitflag = 'fiteach';
end
debug = 0;
if debug
    debugflag = 'yesBust';
else
    debugflag = 'noBust';
end
welch = 1;
if welch
    spectrumflag = 'welch';
else
    spectrumflag = 'periodogram';
end
overlap = 1;
if overlap
    overlapflag = '50overlap';
else
    overlapflag = '0overlap';
end

windowed = 1;

% power or amplitude
factor = 1; %0.5
% figure format
format = 'png';

%creat folders for figures
op_path = char(strjoin({op_path, 'spectral', ...
    normflag, debugflag, overlapflag, fitflag, filesep}, filesep));
mkdir = strjoin({'mkdir -p', op_path});
system(char(mkdir));

% for i_dtau = 1:length(v_dtau)
%     dtau = v_dtau(i_dtau);
%     [tslist, telist] = dividet_v2(t, dtau, 10);
%     [tslist, telist] = dividet_v1(t, dtau, 10);
%     [tslist, telist] = dividet_v3(t, dtau*3/4, 10);
%     [tslist telist telist - tslist]
%     for tt = 1:length(telist)
if ~isnan(vmag) && ~isnan(vang)
    Fs = 100;
    
    %get time-shifted signals
    %     [ph_truncated, pwr_truncated] = ...
    %         plot_shifted(xdata, tpeak, rcvr_op, sitenum_op, combos, Fs);
    for rr = 1:size(rcvr_op, 1)
        if strcmp(rcvr_op(rr, :), 'ASTRArx') && strcmp(doy, '342') && strcmp(year, '2013')
            AZ(tt, rr) = mean(AZ(tt, [1:rr - 1, rr + 1:end]));
            ZE(tt, rr) = mean(ZE(tt, [1:rr - 1, rr + 1:end]));
        end
        zmin = 100e3;
        zmax = 1000e3;
        Lmin = 25e3;
        step = 25e3;
        
        %Amplitude and phase of the receivered signal
        
        %original signals
        pwr_o = xdata{rr}(:, 2); %+ 0.25*randn(l,1);
        ph_o = xdata{rr}(:, 3); %+ 0.25*randn(l,1);
        %time-shifted, aligned and truncated signals
        pwr = pwr_o; %+ 0.25*randn(l,1);
        ph = ph_o; %+ 0.25*randn(l,1);
        
        l = length(pwr);
        NFFT = 2^nextpow2(l);
        l = NFFT;
        
        if windowed
            window_welch = [];
            window_period = [];
        else
            window_welch = ones(l, 1);
            window_period = ones(l, 1);
        end
        
        if overlap
            [Spwr_obs_welch{tt, rr}, ~] = pwelch(factor * log(pwr), window_welch, [], l, Fs, 'psd');
            [Sph_obs_welch{tt, rr}, f] = pwelch(ph, window_welch, [], l, Fs, 'psd');
        else
            [Spwr_obs_welch{tt, rr}, ~] = pwelch(factor * log(pwr), window_welch, 0, l, Fs, 'psd');
            [Sph_obs_welch{tt, rr}, f] = pwelch(ph, window_welch, 0, l, Fs, 'psd');
        end
        [Spwr_obs_period{tt, rr}, ~] = periodogram(factor * log(pwr), window_period, l, Fs, 'psd');
        [Sph_obs_period{tt, rr}, f] = periodogram(ph, window_period, l, Fs, 'psd');
        
        R_obs_welch(:, rr) = Spwr_obs_welch{tt, rr} ./ Sph_obs_welch{tt, rr};
        R_obs_period(:, rr) = Spwr_obs_period{tt, rr} ./ Sph_obs_period{tt, rr};
    end
    
    %Welch or Periodogram?
    if welch
        R_obs = R_obs_welch;
    else
        R_obs = R_obs_period;
    end
    
    [Lgrid, zgrid] = meshgrid(Lmin:step:zmax, zmin:step:zmax);
    for rr = 1:size(rcvr_op, 1)
        ep = NaN(size(Lgrid));
        for i = 1:size(Lgrid, 1)
            for j = 1:size(Lgrid, 2)
                L = Lgrid(i, j);
                z = zgrid(i, j);
                if L < z
                    if fitall == 0
                        [R_rytov, k_par, k_par_index, R_Bust] = Lz(vmag, vang, AZ(tt, rr), ZE(tt, rr), L, z, f);
                        R_obs_c = R_obs(k_par_index, rr);
                    else
                        [R_rytov, k_par, k_par_index, R_Bust] = Lz(vmag, vang, AZ(tt, :), ZE(tt, :), L, z, f);
                        R_obs_c = R_obs(k_par_index, :);
                    end
                    
                    R_rytov_c = R_rytov(k_par_index, :);
                    R_Bust_c = R_Bust(k_par_index, :);
                    
                    % normalized?
                    if normalize == 0
                        sumsquared = (R_obs_c - R_rytov_c).^2;
                    else
                        sumsquared = ((R_obs_c - R_rytov_c) ./ R_obs_c).^2;
                    end
                    
                    epsqr = mean(sumsquared(:), 'omitnan');
                    ep(i, j) = epsqr;
                    %                             if epsqr > 10
                    %                                 ep(i,j) = NaN;
                    %                                 continue;
                    %                             end
                end
            end
        end
        figj = figure;
        
        ep_min(:, rr) = min(min(ep));
        ep_max(:, rr) = max(max(ep));
        L_hat(:, rr) = Lgrid(ep == ep_min(:, rr));
        z_hat(:, rr) = zgrid(ep == ep_min(:, rr));
        %                 mesh(Lgrid/10^3, zgrid/10^3, ep, ep);
        pcolor(Lgrid/10^3, zgrid/10^3, ep);
        %                 shading(gca, 'flat');
        %                 set(gca, 'layer', 'top');
        caxis([0, 5]);
        hold on;
        %                 plot3(L_hat/10^3,z_hat/10^3,ep_min,'ro');
        plot(L_hat(:, rr)/10^3, z_hat(:, rr)/10^3, 'ro');
        xlabel('Thickness $L$ [km] ');
        ylabel('Top height $z$ [km]');
        zlabel('$\epsilon^2$');
        title(['$\hat{L} = $', num2str(L_hat(:, rr)/10^3), ', ', ...
            '$\hat{z} = $', num2str(z_hat(:, rr)/10^3), ', ', ...
            '$\epsilon^2_{min} = $', num2str(ep_min(:, rr)), ', ', ...
            'for ', sitenum_op{rr, :}]);
        %                 set(gca, 'Zscale', 'log');
        %         zlim([0 10]);
        %                 view([45,15]);
        tightfig;
        cb = colorbar;
        set(get(cb, 'YLabel'), 'String', '$\epsilon^2$', 'interpreter', 'latex');
        plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_CostFunction_', ...
            num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
            datestr(init_time, 'HHMM'), 'UT_', num2str(factor), '_', num2str(zmax)];
        plotpath = [op_path, plotname];
        saveas(gcf, plotpath, format);
        close;
    end
    xl = [-inf, inf]; %xl = [1e-3 1e-1];
    %             return;
    for rr = 1:size(rcvr_op, 1)
        figall = figure;
        if welch
            loglog(k_par, Spwr_obs_welch{tt, rr}, 'b', ...
                k_par, Sph_obs_welch{tt, rr}, 'k', ...
                k_par, R_obs_welch(:, rr), 'r');
        else
            %                     hold on;
            loglog(k_par, Spwr_obs_period{tt, rr}, 'c', ...
                k_par, Sph_obs_period{tt, rr}, 'g', ...
                k_par, R_obs_period(:, rr), 'm');
        end
        xlim(xl);
        %                     set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
        title('Observed Log-Amplitude to Phase Power Spectrum Ratio');
        legend({'Log-Amplitude', 'Phase', 'Ratio'}, 'location', 'southwest');
        xlabel(['Wavenumber along Drift Velocity Direction $k_{\parallel}$ [rad/m], ', sitenum_op{rr, :}]);
        %                     legend({'Phase','Log_{10} Power'},'location','best')
        tightfig;
        plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_ObservedRatio_', ...
            num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
            datestr(init_time, 'HHMM'), 'UT'];
        plotpath = [op_path, plotname];
        saveas(gcf, plotpath, format);
        close;
    end
    fighat = figure;
    if size(rcvr_op, 1) > 2
        set(gcf, 'papersize', [8, 2 * size(rcvr_op, 1)], ...
            'paperposition', [0, 0, 8, 2 * size(rcvr_op, 1)], ...
            'paperpositionmode', 'auto', ...
            'position', [0, 0, 8, 2 * size(rcvr_op, 1)]);
    end
    [sp, ~] = tight_subplot(size(rcvr_op, 1), 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
    for rr = 1:size(rcvr_op, 1)
        [R_rytov_hat, k_par, k_par_index, R_Bust] = Lz(vmag, vang, ...
            AZ(tt, rr), ZE(tt, rr), L_hat(:, rr), z_hat(:, rr), f);
        [R_rytov_fixed, ~, k_par_index_Bust, R_Bust_fixed] = Lz(vmag, vang, ...
            AZ(tt, rr), ZE(tt, rr), 500e3, 200e3, f);
        k_par_c = k_par(k_par_index, :);
        R_rytov_hat_c(:, rr) = R_rytov_hat(k_par_index, :);
        R_obs_c(:, rr) = R_obs(k_par_index, rr);
        if debug
            loglog(sp(rr), k_par_c, R_obs_c(:, rr), 'r', ...
                k_par_c, R_rytov_hat_c(:, rr), 'c', ...
                k_par_c, R_Bust(k_par_index), 'k', ...
                k_par_c, R_rytov_fixed(k_par_index), 'g', ...
                k_par_c, R_Bust_fixed(k_par_index), 'b');
            legend(sp(rr), ['Observed, ', sitenum_op{rr, :}], ...
                ['R, $\hat{L}=', num2str(L_hat(:, rr)/10^3), '$, $\hat{z}=', num2str(z_hat(:, rr)/10^3), '$'], ...
                'B, -', ...
                'R, $\hat{L}=500$, $\hat{z}=200$', ...
                'B, -', ...
                'location', 'northwest', 'orientation', 'horizontal');
        else
            loglog(sp(rr), k_par_c, R_obs_c(:, rr), 'r', ...
                k_par_c, R_rytov_hat_c(:, rr), 'c');
            legend(sp(rr), ['Observed, ', sitenum_op{rr, :}], ...
                ['Rytov, $\hat{L}=', num2str(L_hat(:, rr)/10^3), ...
                '$, $\hat{z}=', num2str(z_hat(:, rr)/10^3), '$'], ...
                'location', 'southeast');
        end
        xlim(sp(rr), xl);
        ylim(sp(rr), [10^-5, 10^2.5]);
        %                     set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
        %                     ylim([1e-5*0.99 1e2*1.01]);
        MEGA_LZ(tt, rr, :) = [datenum(tslist(tt)/24/3600+init_time), ...
            datenum(telist(tt)/24/3600+init_time), ...
            L_hat(:, rr) / 10^3, z_hat(:, rr) / 10^3, ep_min(:, rr)];
    end
    MEGA_LZ
    title(sp(1), 'Rytov and Observed Log-Amplitude to Phase Power Spectrum Ratio');
    xlabel(sp(rr), 'Wavenumber along Drift Velocity Direction $k_{\parallel}$ [rad/m]');
    set(sp(1:(rr - 1)), 'xticklabel', []);
    set(sp(2:2:rr), 'yticklabel', []);
    tightfig;
    plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_RytovObserved_', ...
        num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
        datestr(init_time, 'HHMM'), 'UT_', num2str(factor), '_', num2str(zmax)];
    plotpath = [op_path, plotname];
    saveas(gcf, plotpath, format);
    %             print(gcf, [plotpath, '_pver.png'], '-dpng', '-r600');
    %                 close;
    
    close;
else
    for rr = 1:size(rcvr_op, 1)
        MEGA_LZ(tt, rr, :) = [datenum(tslist(tt)/24/3600+init_time), ...
            datenum(telist(tt)/24/3600+init_time), ...
            NaN, NaN, NaN];
    end
end

%     end
%     if isempty(tslist) && isempty(telist)
%         for rr = 1:size(rcvr_op, 1)
%             MEGA_LZ(tt,rr,:) = [NaN; NaN; NaN; NaN; NaN];
%         end
%     end
% end
