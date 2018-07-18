% Script to test effect of noise on estimates of tau_a and tau_c.
%
% Seebany Datta-Barua
% 22 Jan 2017

clear
close all
tic
% Your variables are:
%
% init_t_utc  rcvr_op     sitenum_op  xdata_PRN
% prn         signal      tlim
% prn = 23
% sitenum_op =
%
%     'IIT-1'
%     'IIT-3'
%     'IIT-11'
%     'IIT-13'
%     'IIT-15'
load('../HRdataforDrBust/Bust_PRN23_2013_342');

for i = 1:size(sitenum_op)
    rxdata = xdata_PRN{i};
    t(:, i) = rxdata(:, 1);
    phi(:, i) = rxdata(:, 3);
end

% Generate different ensemble members.
% Number of ensemble members in the simulation.
num_ensemble = 100;
noise = randn(numel(phi(:, 1)), numel(sitenum_op), num_ensemble);
t0 = min(t(1, :));

% figure(1)
% subplot(211)
% % First plot data without noise
% plot(t-t0,phi)
% legend(sitenum_op)
% title('Unaltered phase data from SAGA to PRN 23')
% ylabel('Detrended filtered 100 Hz phase [rad]')
%
% % Second plot is data with noise.
% subplot(212)
% plot(t-t0,phi+noise(:,:,1))
% title('One ensemble member of SAGA data with added noise')
% xlabel(['Seconds since 03:' num2str(2620/60) ' UT'])
% ylabel('Detrended filtered 100 Hz phase [rad]')
%

for k = 1:num_ensemble
    % Compute auto- and cross-correlations.
    rho = xcorr(phi);
    rhon = xcorr(phi+noise(:, :, k));
    
    % The columns of rho that contain auto-correlations.
    num_rxs = numel(sitenum_op);
    acorrcols = ([1:num_rxs] - 1) * num_rxs + ([1:num_rxs] - 1) + 1;
    maxrhoii = max(rho(:, acorrcols));
    % Noisy version
    maxrhoiin = max(rhon(:, acorrcols));
    
    % Initialize sample observations
    y = [];
    yn = [];
    diffy = [];
    
    for i = 1:1 %numel(sitenum_op)-1
        tau = [-(t(end:-1:1, i) - t0); t(2:end, i) - t0];
        
        %     for j = acorrcols(i)+1:acorrcols(i+1)
        for j = i + 1:numel(sitenum_op)
            %         j
            normfactor = sqrt(maxrhoii(i)*maxrhoii(j));
            % Noisy version
            normfactorn = sqrt(maxrhoiin(i)*maxrhoiin(j));
            
            % Normalize the autocorr and xcorr functions.
            normacorr = rho(:, acorrcols(i)) / normfactor;
            normxcorr = rho(:, acorrcols(i)+j-1) / normfactor;
            % Noisy.  Not sure it's right to normalize by the noise-free
            % version but the noisy version has a spike of almost double the
            % true peak value.
            normacorrn = rhon(:, acorrcols(i)) / normfactor;
            normxcorrn = rhon(:, acorrcols(i)+j-1) / normfactor;
            
            
            % Pick one side of the main lobe above rho_cutoff.
            rho_cutoff = 0.6;
            arows = find(normacorr > rho_cutoff & tau >= 0);
            xrows = find(normxcorr > rho_cutoff & ...
                tau >= tau(find(normxcorr == max(normxcorr))));
            
            shortacorr = normacorr(arows);
            shortatau = tau(arows);
            shortxcorr = normxcorr(xrows);
            shortxtau = tau(xrows);
            
            %%% Noisy version
            % Pick one side of the main lobe above rho_cutoff.
            rho_cutoffn = rho_cutoff;
            arowsn = find(normacorrn > rho_cutoffn & tau >= 0);
            xrowsn = find(normxcorrn > rho_cutoffn & ...
                tau >= tau(find(normxcorrn == max(normxcorrn))));
            
            shortacorrn = normacorrn(arowsn);
            shortataun = tau(arowsn);
            shortxcorrn = normxcorrn(xrowsn);
            shortxtaun = tau(xrowsn);
            
            figure(2)
            % First plot is one set of autocorrelations, noisefree and  noisy.
            subplot(211)
            plot(tau, normacorrn, 'g', tau, normacorr, 'b') %, ...
            %             tau(arows), shortacorr, 'r', tau(arowsn), shortacorrn,'c')
            legend('noise-added', 'raw data') %,'raw, above \rho_c', ...
            title(['Auto-correlations']) % of data from ' sitenum_op{i}])
            %         legend('raw data','noise-added', 'raw, above \rho_c', ...
            %             'noisy, above \rho_c')
            %         title(['Auto-correlation of data from ' sitenum_op{i}])
            %         ax = axis;
            %         axis([-4 4 ax(3) ax(4)])
            hold on
            
            % Second plot is one pair of cross-correlations, noisefree & noisy.
            figure(2)
            subplot(212)
            plot(tau, normxcorrn, 'g', tau, normxcorr, 'b') %,...
            %             tau(xrows), shortxcorr, 'r', tau(xrowsn), shortxcorrn,'c')
            legend('noise-added', 'raw data') %,'raw, above \rho_c', ...
            title(['Cross-correlations']) % of data from ' sitenum_op{i}])
            %         legend('raw data','noise-added', 'raw, above \rho_c', ...
            %             'noisy, above \rho_c')
            %         title(['Cross-correlation of data from ' sitenum_op{i} ' and ' ...
            %             sitenum_op{j}])
            %         ax = axis;
            %         axis([-4 4 ax(3) ax(4)])
            hold on
            
            % Repeat matrix of each value of cross-correlation.
            xcorrmat = repmat(shortxcorr', size(shortacorr));
            
            % Difference from the autocorrelations.
            acorrmat = repmat(shortacorr, size(shortxcorr'));
            absdiffmat = abs(xcorrmat-acorrmat);
            [minrows, mincols] = find(absdiffmat == repmat(min(absdiffmat), ...
                size(absdiffmat, 1), 1));
            tau_a = shortatau(minrows);
            tau_x = shortxtau;
            
            % Repeat matrix of each value of cross-correlation.
            xcorrmatn = repmat(shortxcorrn', size(shortacorrn));
            
            % Difference from the autocorrelations.
            acorrmatn = repmat(shortacorrn, size(shortxcorrn'));
            absdiffmatn = abs(xcorrmatn-acorrmatn);
            [minrowsn, mincolsn] = find(absdiffmatn == repmat(min(absdiffmatn), ...
                size(absdiffmatn, 1), 1));
            tau_an = shortataun(minrowsn);
            tau_xn = shortxtaun;
            
            % Observations are differences of the squares.
            % Truncate differences.
            len = min(numel(tau_an), numel(tau_a));
            diffy = [diffy; tau_an(1:len).^2 - tau_xn(1:len).^2, ...
                - (tau_a(1:len).^2 - tau_x(1:len).^2)];
            y = [y; tau_a.^2 - tau_x.^2];
            yn = [yn; tau_an(1:len).^2 - tau_xn(1:len).^2];
            figure(4)
            plot(diffy)
            hold on
        end
    end
    
    figure(3)
    hold on
    plot(yn, 'g')
    title('Observation array')
    %         legend('noise-free','noisy')
    ensemble_yn{k} = yn;
end

%Complete Fig 3 with the noise-free y.
figure(3)
plot(y)
xlabel('Element number')
ylabel('\tau_a^2 - \tau_c^2 [s]')

figure(2)
subplot(212)
ax = axis;
axis([-4, 4, ax(3), ax(4)])
xlabel('Lag \tau  [s]')
grid on
ylabel('Normalized auto-correlation')

subplot(211)
grid on
ax = axis;
axis([-4, 4, ax(3), ax(4)])
ylabel('Normalized cross-correlation')
toc