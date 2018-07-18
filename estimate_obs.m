function [tauaarrn_, taucarrn_, ccvalarrn_, ccerrarrn_, ...
    tauaarr_, taucarr_, ccvalarr_, ccerrarr_] = ...
    estimate_obs(rcvr_op, xdata, combos, flag)
% Construct observations for the linear system
dbstop if error;

% rng('default');
tic;
if nargin == 0
    disp('debug begins');
    close all;
    clear;
    load('testdata1.mat');
    %     load('../HRdataforDrBust/Bust_PRN23_2013_342');
    %     xdata = xdata_PRN;
else
    save('testdata.mat');
    %     exit;
end

% Specify which phase data are used
if strcmp(flag, 'phase')
    colnum = 3;
    dt = 0.01;
elseif strcmp(flag, 'power')
    colnum = 2;
    dt = 0.01;
elseif strcmp(flag, 'JROph')
    colnum = 3;
    dt = 0.02;
elseif strcmp(flag, 'JROpow')
    colnum = 2;
    dt = 0.02;
end

data = horzcat(xdata{:});
tmat = data(:, 1:size(xdata{1}, 2):end);
obsmat = data(:, colnum:size(xdata{1}, 2):end);

[rhoarr, lag] = xcorr(obsmat);
rhomax = max(rhoarr);

% Permutations of receiver pairs including auto-correlation
perms = sortrows([combos; fliplr(combos); repmat([1:size(rcvr_op, 1)]', 1, 2)]);
% [~, accols] = ismember(repmat([1:size(rcvr_op, 1)]', 1, 2), perms, 'rows');

% Generate normally distributed gaussian noise
numsim = 1000;
scale = mean(std(obsmat, 0, 1));
scale = 0.25;
noisearr = scale * randn(size(obsmat, 1), size(obsmat, 2), numsim);
if nargin == 0
    figphase = figure;
    subplot(2, 1, 2);
    hold on;
    hwn = plot(tmat, obsmat+noisearr(:, :, 1), 'g.');
    hwon = plot(tmat, obsmat, 'k', 'linewidth', 0.5);
    legend([hwn(1), hwon(1)], {'noisy $\tilde{\phi}$', 'original $\phi$'}, ...
        'orientation', 'horizontal', 'location', 'northwest');
    ylim([-2 * pi, 2 * pi]);
    xlabel('Time [s]');
    ylabel('Phase [rad]');
    xlim([min(min(tmat)) - 1, max(max(tmat)) + 1]);
    title({'Filtered receiver phase measurements ensemble n = 1'; ['with white gaussian noise $\sigma_w^2 = $', num2str(scale^2)]});
    [~, path] = ver_chk;
    saveas(figphase, [path, 'IQ_zoom_', ...
        num2str(min(tmat(:, 1)), '%.f'), '_', num2str(max(tmat(:, 1)), '%.f'), 'a.eps'], 'epsc2');
    close(figphase);
end

[taucarrn, tauaarrn, ccvalarrn, ccerrarrn] = deal(cell(length(combos), numsim));
[taucarr, tauaarr, ccvalarr, ccerrarr] = deal(cell(length(combos), 1));
figobs = figure;
for k = 1:numsim
    [tobsn_k, tobs_k] = deal([]);
    %     if verLessThan('matlab', '9.1.0')
    %         noisearr = randn(size(obsmat));
    %     else
    %         noisearr = wgn(size(obsmat, 1), size(obsmat, 2), 0);
    %     end
    rhoarrn = xcorr(obsmat+noisearr(:, :, k));
    rhomaxn = max(rhoarrn);
    %loop through pairs of receivers
    for i = 1:length(combos)
        [~, ijcol] = ismember(combos(i, :), perms, 'rows');
        %         [~, jicol] = ismember(combos(i,:), perms, 'rows');
        [~, iicol] = ismember([combos(i, 1), combos(i, 1)], perms, 'rows');
        [~, jjcol] = ismember([combos(i, 2), combos(i, 2)], perms, 'rows');
        
        %Normalization
        normscale = sqrt(rhomax(iicol)*rhomax(jjcol));
        %         abc = rhoarr(:,iicol) + ...
        %             xcorr(noisearr(:,combos(i, 1),k), noisearr(:,combos(i, 1),k)) + ...
        %             xcorr(noisearr(:,combos(i, 1),k), obsmat(:,combos(i, 1))) + ...
        %             xcorr(obsmat(:,combos(i, 1)), noisearr(:,combos(i, 1),k));
        normscale1 = rhomax(iicol);
        normscale1n = rhomaxn(iicol);
        normscale2 = rhomax(ijcol);
        normscale2n = rhomaxn(ijcol);
        normscalen = sqrt(rhomaxn(iicol)*rhomaxn(jjcol));
        cc = rhoarr(:, ijcol) / normscale;
        %     cc2(:, i) = rhomat(:,jicol) / normscale;
        ac = rhoarr(:, iicol) / normscale;
        %     ac2(:, i) = rhomat(:,jjcol) / normscale;
        
        ccn = rhoarrn(:, ijcol) / normscale;
        %     cc2n(:, i) = rhonmat(:,jicol) / normscale;
        acn = rhoarrn(:, iicol) / normscale;
        %         acn(lag==0,i) = NaN;
        %         ac(lag==0,i) = NaN;
        %     ac2n(:, i) = rhonmat(:,jjcol) / normscale;
        
        %     B = 1000e-4;
        %     rhoiin = rhomaxn(iicol) * sinc(2*pi*B*lag'*dt);
        %     plot(lag*dt, rhon(:, iicol), lag*dt, rhoiin)
        
        %     %average auto-correlation of pairwise receiver signals
        %     ac(:, i) = (ac(:, i) + ac2(:, i)) / 2;
        
        %     hold on;
        
        
        if nargin == 0 && i == 1 && k == 1
            figrawcorr = figure;
            hold on;
            h1 = plot(lag*dt, cc, 'k', lag*dt, ac, 'c');
            [ccl, ccr] = findmainlobe(cc, lag);
            [ccnl, ccnr] = findmainlobe(ccn, lag);
            title(['Receiver ', num2str(combos(i, 1)), ' \& ', num2str(combos(i, 2))]);
            h2 = plot(lag*dt, ccn, 'g.', lag*dt, acn, 'r.', 'markersize', 1);
            xlabel('Lag [s]');
            [acl, acr] = findmainlobe(ac, lag);
            [acnl, acnr] = findmainlobe(acn, lag);
            legend('$\rho_{ij}$', '$\rho_{ii}$', ...
                '$\tilde{\rho_{ij}}$', '$\tilde{\rho_{ii}}$');
            xlim([min([acl, ccl, acnl, ccnl]), max([acr, ccr, acnr, ccnr])]*dt);
            ylim([-0.08, 1.02]);
            tightfig;
            saveas(figrawcorr, '../ccexp.pdf');
            %             keyboard;
            close(figrawcorr);
        end
        
        if nargin == 0 && k == 1 && i == 1
            debugflag = 0;
        else
            debugflag = 1;
        end
        if k == 1
            [rowscc{i}, tobs, tau_c, tau_a, ccval, ccerr] = ...
                findrows(cc, ac, lag, dt, 'original', [], debugflag);
        end
        
        [~, tobsn, tau_cn, tau_an, ccvaln, ccerrn] = ...
            findrows(ccn, acn, lag, dt, 'noisy', rowscc{i}, debugflag);
        
        tobs_k = [tobs_k, tobs];
        tobsn_k = [tobsn_k, tobsn];
        
        taucarrn{i, k} = tau_cn;
        tauaarrn{i, k} = tau_an;
        ccvalarrn{i, k} = ccvaln;
        ccerrarrn{i, k} = ccerrn;
        
        if k == 1
            taucarr{i} = tau_c;
            tauaarr{i} = tau_a;
            ccvalarr{i} = ccval;
            ccerrarr{i} = ccerr;
        end
    end
    
    if nargin == 1
        figure(figobs);
        hold on;
        if k == 1
            tobsarr = zeros(numsim, length(tobs_k));
            errarr = zeros(numsim, length(tobs_k));
            plot(tobs_k, 'k.-', 'linewidth', 2);
        end
        errarr(k, :) = (tobsn_k - tobs_k)';
        plot(tobsarr(k, :), 'g');
        plot(errarr(k, :), 'r.');
        drawnow;
    end
end
% close(figobs);
[taucarrn_, tauaarrn_, ccvalarrn_, ccerrarrn_, ...
    taucarr_, tauaarr_, ccvalarr_, ccerrarr_] = deal(cell(length(combos), 1));
for i = 1:length(combos)
    taucarrn_{i, :} = vertcat(taucarrn{i, :})';
    tauaarrn_{i, :} = vertcat(tauaarrn{i, :})';
    ccvalarrn_{i, :} = vertcat(ccvalarrn{i, :})';
    ccerrarrn_{i, :} = vertcat(ccerrarrn{i, :})';
    taucarr_{i, :} = horzcat(taucarr{i, :})';
    tauaarr_{i, :} = horzcat(tauaarr{i, :})';
    ccvalarr_{i, :} = horzcat(ccvalarr{i, :})';
    ccerrarr_{i, :} = horzcat(ccerrarr{i, :})';
end
toc;
end

function [tmprowscc, tobs, tau_c, tau_a, ccval, errs] = ...
    findrows(cc, ac, lag, dt, flag, rowscc, debugflag)
% figpks = figure;
% findpeaks(cc, lag, 'annotate', 'extents', 'widthreference', 'halfheight');
% [pks, loc, ~, ~] = ...
%     findpeaks(cc, lag, 'annotate', 'extents', 'widthreference', 'halfheight');

%above a cut-off, i.e. 0.6
rhocutoff = 0.5;
nearzero = 1e-2;
% nearzero = 0.1;

lagmax = lag(cc == max(cc));
[laglcc, lagrcc] = findmainlobe(cc, lag);
[laglac, lagrac] = findmainlobe(ac, lag);
tmprowsac = find(ac' > 0 & lag > 0 & lag < lagrac);

if ~strcmp(flag, 'noisy')
    %take the lag positive side of the main lobe
    tmprowscc = find(cc' >= rhocutoff & lag >= lagmax & lag < lagrcc);
else
    tmprowscc = rowscc;
end

if debugflag == 0
    figdebug = figure;
    %     ac(lag==0) = NaN;
    %     plot(lag*dt, cc, 'k', lag*dt, ac, 'k');
    hold on;
    %     plot(laglac*dt,ac(lag==laglac),'o');
    %     plot(lagrac*dt,ac(lag==lagrac),'o');
    hc = plot(lag*dt, cc, 'r');
    %     plot(lag(tmprowscc_)*dt, cc(tmprowscc_), 'rs-');
    ha = plot(lag(tmprowsac)*dt, ac(tmprowsac), 'c');
    if strcmp(flag, 'noisy')
        rhostr = '\tilde{\rho}';
    else
        rhostr = '\rho';
    end
    title([flag, ', $\rho_{cutoff}$ = ', num2str(rhocutoff), ...
        ', $\min |', rhostr, '_{ij}-', rhostr, '_{ii}| \le$', num2str(nearzero)]);
    legend([hc, ha], {['$', rhostr, '_{ij}$'], ['$', rhostr, '_{ii}$']});
    xlabel('Lag [s]');
    if ~isempty(tmprowscc)
        xlim(dt*[min(lag(tmprowscc(1)), 0), max(lagrcc, lagrac)]);
    else
        xlim(dt*[min(lagmax, 0), max(lagrcc, lagrac)]);
    end
    %     keyboard;
    %         close;
end

if ~isempty(tmprowscc)
    
    ccmat = repmat(cc(tmprowscc)', size(ac(tmprowsac)));
    acmat = repmat(ac(tmprowsac), size(cc(tmprowscc)'));
    [errs, rowsmin] = min(abs(ccmat-acmat), [], 1);
    
    %remove any cross-auto pair with large errors
    invalid = errs >= nearzero;
    %     rowsmin(errs >= nearzero) = NaN;
    %     tmprowscc(errs >= nearzero) = NaN;
    %     errs(errs >= nearzero) = NaN;
    
    
    %make sure one tau_a corresponds to one tau_c
    if debugflag == 0
        for indcc = find(~isnan(tmprowscc) & ~invalid)
            plot([lag(tmprowscc(indcc)) * dt, lag(tmprowsac(rowsmin(indcc))) * dt], ...
                [cc(tmprowscc(indcc)), ac(tmprowsac(rowsmin(indcc)))], 'k.-', ...
                'linewidth', 0.5);
        end
        %         for rowac = unique(rowsmin, 'stable')
        %             rowunique = find(rowsmin == rowac);
        %             h1 = plot(lag(tmprowscc(rowunique))*dt, cc(tmprowscc(rowunique)), 'k.');
        %             rowminerr = rowunique(errs(rowunique) == min(errs(rowunique)));
        %             h2 = plot(lag(tmprowscc(rowminerr))*dt, cc(tmprowscc(rowminerr)), 'k*');
        %             h3 = plot([lag(tmprowscc(rowminerr))*dt lag(tmprowsac(rowac))*dt], ...
        %                 [cc(tmprowscc(rowminerr)) ac(tmprowsac(rowac))], 'k.-');
        %             errs(rowunique(rowunique ~= rowminerr)) = [];
        %             rowsmin(rowunique(rowunique ~= rowminerr)) = [];
        %             tmprowscc(rowunique(rowunique ~= rowminerr)) = [];
        %         end
        tightfig;
        saveas(figdebug, ['../getobs', flag, '.pdf']);
        %         keyboard;
        close(figdebug);
    end
    
    if any(invalid)
        %         disp('There are NaNs in $tau_a$');
    end
    tau_a = lag(tmprowsac(rowsmin)) * dt;
    tau_a(invalid) = NaN;
    if isempty(tau_a)
        tau_a = NaN(1, length(tmprowscc));
        errs = tau_a;
    end
    tau_c = lag(tmprowscc) * dt;
    ccval = cc(tmprowscc)';
    tobs = tau_a.^2 - tau_c.^2;
    if any(abs(tau_a-tau_c) > 10)
        lagmax
        [laglcc, lagrcc] = findmainlobe(cc, lag)
        [laglac, lagrac] = findmainlobe(ac, lag)
    end
    if debugflag == 0
        figcc = figure;
        hold on;
        plot(lag*dt, cc, lag*dt, ac);
        plot(lag(tmprowscc)*dt, cc(tmprowscc), 'r', 'LineWidth', 2);
        plot(lag(tmprowsac)*dt, ac(tmprowsac), 'c', 'LineWidth', 2);
        plot(tau_c, cc(tmprowscc), 'k.');
        plot(tau_a, ac(tmprowsac(rowsmin)), 'k.');
        close(figcc);
    end
else
    tobs = [];
    tau_c = [];
    tau_a = [];
    ccval = [];
    errs = [];
end
end

function [lagmainl, lagmainr] = findmainlobe(cc, lag)
if all(abs(flip(cc)-cc) < 1e-13)
    lagmax = 0;
else
    lagmax = lag(cc == max(cc));
end

lagzerosl = lag((circshift(cc, 1) .* cc < 0 | circshift(cc, -1) .* cc < 0) & lag' < lagmax);
lagzerosr = lag((circshift(cc, 1) .* cc < 0 | circshift(cc, -1) .* cc < 0) & lag' > lagmax);

[~, indminr] = min(lagzerosr-lagmax);
lagmainr = lagzerosr(indminr);
[~, indminl] = max(lagzerosl-lagmax);
lagmainl = lagzerosl(indminl);
end
