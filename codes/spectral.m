%Spectral Analysis Part
load lz.mat
v_ccmin = [0.6];
v_dtau = 60;
for i_dtau = 1:length(v_dtau)
    dtau = v_dtau(i_dtau);
    [tslist, telist] = dividet_v2(t, dtau, 10);
    [tslist, telist] = dividet_v1(t, dtau, 10);
%     [tslist, telist] = dividet_v3(t, dtau*3/4, 10);
    [tslist telist telist - tslist]
    for tt = 1:length(telist)        
        if ~isnan(vmag) && ~isnan(vang) && 0
            %spectral analysis
            %                     if (strcmp(doy,'051') || strcmp(doy,'342')) && ccmin == 0.6 && dtau == 60
            Fs = 100;
            for rr = 1:size(rcvr_op, 1)
                if strcmp(rcvr_op(rr,:), 'ASTRArx') && strcmp(doy, '342') && strcmp(year, '2013')
                    AZ(tt, rr) = mean(AZ(tt, [1:rr - 1, rr + 1:end]));
                    ZE(tt, rr) = mean(ZE(tt, [1:rr - 1, rr + 1:end]));
                end
                zmin = 150e3; zmax = 500e3; Lmin = 0e3; step = 10e3;
                %Amplitude and phase of the receivered signal
                l = length(xdata{rr}(:, 1));
                pwr = xdata{rr}(:, 2); %+ 0.25*randn(l,1);
                ph = xdata{rr}(:, 3); %+ 0.25*randn(l,1);        
                NFFT = 2^nextpow2(l);
               [Spwr_obs{tt, rr}, ~] = pwelch(log(pwr), [], [], NFFT, 'power');
        %         [Spwr_obs{tt, rr}, f] = periodogram(pwr, [], [], NFFT, 'psd');
                [Sph_obs{tt, rr}, ~] = pwelch(ph, [], [], NFFT, 'power');
        %         [Spwr_obs{tt, rr}, f] = periodogram(pwr, [], [], NFFT, 'psd');
                R_obs = Spwr_obs{tt, rr} ./ Sph_obs{tt, rr};
                [Lgrid, zgrid] = meshgrid(Lmin:step:zmax,zmin:step:zmax);
                ep = NaN(size(Lgrid));
                for i = 1:size(Lgrid,1)
                    for j = 1:size(Lgrid,2)
                        L = Lgrid(i,j);
                        z = zgrid(i,j);
                        if L < z
                            [R_rytov, ~, k_par_index] = Lz(Fs, NFFT, vmag, vang, AZ(tt, rr), ZE(tt, rr), L/10^3, z/10^3);
                            R_rytov_c = R_rytov(k_par_index,:);
                            R_obs_c = R_obs(k_par_index,:);
                            epsqr = mean(((R_rytov_c - R_obs_c) ./ R_obs_c) .^2,'omitnan');                     
                            ep(i,j) = sqrt(epsqr);
                        else
                            ep(i,j) = NaN;
                        end
                    end
                end

                figj = figure;

                ep_min = min(min(ep));
                L_hat = Lgrid(ep==ep_min);
                z_hat = zgrid(ep==ep_min);

                mesh(Lgrid/10^3, zgrid/10^3, ep, ep);
                hold on;
                plot3(L_hat/10^3,z_hat/10^3,ep_min,'ro');
                xlabel('Thickness L[km] '); ylabel('Top height z [km]'); zlabel('$\epsilon$');
                title(['$\hat{L} =$', num2str(L_hat/10^3), ', ', ...
                    '$\hat{z} =$', num2str(z_hat/10^3), ', ', ...
                    '$\epsilon_{min} =$', num2str(ep_min)]);
                set(gca, 'Zscale', 'log');
        %         zlim([0 10]);
                view([45,15]);
                plotname = ['PRN', num2str(prn), '_', sitenum_op{rr,:}, '_CostFunction_', ...
                    num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                    datestr(init_time, 'HHMM'), 'UT'];
                plotpath = [op_path, plotname, '.eps'];
                saveas(gcf, plotpath, 'epsc2');
                close;

                figall = figure;
                [R_rytov, k_par, k_par_index] = Lz(Fs, NFFT, vmag, vang, AZ(tt, rr), ZE(tt, rr), L_hat, z_hat);
                k_par_c = k_par(k_par_index,:);
                R_rytov_c = R_rytov(k_par_index,:);
                R_obs_c = R_obs(k_par_index,:);

                xl = [-inf, inf]; %xl = [1e-3 1e-1];
                loglog(k_par, Spwr_obs{tt, rr}, 'b', k_par, Sph_obs{tt, rr}, 'k', k_par, R_obs, 'r');
                xlim(xl);
                %                     set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
                title('Observed Log-Amplitude to Phase Power Spectrum Ratio');
                legend({'Log-Amplitude', 'Phase', 'Ratio'}, 'location', 'northeast');
                xlabel(['Wavenumber along Drift Velocity Direction $k_{\parallel}$ [rad/m], ', sitenum_op{rr,:}]);
                %                     legend({'Phase','Log_{10} Power'},'location','best')
                plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr,:}, '_ObservedRatio_', ...
                    num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                    datestr(init_time, 'HHMM'), 'UT'];
                plotpath = [op_path, plotname, '.eps'];
                saveas(gcf, plotpath, 'epsc2');
                close;

                fighat = figure;
                loglog(k_par_c, R_obs_c, 'r', k_par_c, R_rytov_c, 'c');
                title('Rytov and Observed Log-Amplitude to Phase Power Spectrum Ratio');
                legend(['Observed,', sitenum_op{rr,:}], ...
                    ['Rytov, $\hat{L} =$ ', num2str(L_hat/10^3), 'km, $\hat{z} =$ ', num2str(z_hat/10^3), 'km'] ...
                    , 'location', 'southeast');
                xlabel('Wavenumber along Drift Velocity Direction $k_{\parallel}$ [rad/m]');
                xlim(xl);
                %                     set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
                %                     ylim([1e-5*0.99 1e2*1.01]);
                plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr,:}, '_RytovObserved_', ...
                    num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                    datestr(init_time, 'HHMM'), 'UT'];
                plotpath = [op_path, plotname, '.eps'];
                saveas(gcf, plotpath, 'epsc2');
                close;
                MEGA_LZ(tt,rr,:) = [datenum(tslist(tt)/24/3600+init_time), ...
                    datenum(telist(tt)/24/3600+init_time),...
                    L_hat / 10^3, z_hat / 10^3, ep_min];
            end
            MEGA_LZ
        else
            for rr = 1:size(rcvr_op, 1)
                MEGA_LZ(tt,rr,:) = [datenum(tslist(tt)/24/3600+init_time), ...
                    datenum(telist(tt)/24/3600+init_time),...
                    NaN, NaN, NaN];
            end    
        end
        
%         return;
    end
    if isempty(tslist) && isempty(telist)
        for rr = 1:size(rcvr_op, 1)
            MEGA_LZ(tt,rr,:) = [NaN; NaN; NaN; NaN; NaN; NaN];
        end
    end
end