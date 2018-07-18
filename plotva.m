function [vmaghat, vanghat, vgehat, vgnhat, arhat, psiahat, vctovhat, argshat] ...
    = plotva(H, Y, CCVAL, CCERR, COMBOS, peak, denu)
%%
if nargin == 0
    clear;
    close all;
    load va1.mat;
    format short g
    dbstop if error;
end
dt = 10e-3;
HYCC = [H, Y, CCVAL];
npt = size(HYCC, 1);
%hardcoded correlation cut off
rho_c = 0.5;
% %hardcoded correlation cut off for SAGA
% rho_c = 0.6;
rhoc_old = 1;
i = 1;

% [ccmax, imax] = max(peak);
% ccmax
% xij = denu(1,imax); yij = denu(2,imax);
% [ccmax2, ii] = max(CCVAL);
% ccmax2
% COMBOS(ii,:)
if max(CCVAL) > rho_c
    while rho_c + (i - 1) * dt < max(CCVAL)
        HYCC_tmp = HYCC(HYCC(:, 7) >= rho_c+i*dt, :);
        ERR_tmp = CCERR(HYCC(:, 7) >= rho_c+i*dt, :);
        rhoc(i, :) = rho_c + i * dt;
        H = HYCC_tmp(:, 1:end-2);
        Y = HYCC_tmp(:, end-1);
        CC = HYCC_tmp(:, end);
        X = pinv(H) * Y;
        if ~isempty(X)
            X = num2cell(X);
            [a, h, b, f, g] = X{:};
            if a * b - h^2 > 0 && a > 0 && b > 0
                [vmag(i, :), vang(i, :), ar(i, :), Psi_a(i, :), Vc(i, :)] = drvel(a, h, b, f, g);
                args(i, :) = [a, h, b, f, g];
                
                %                 [ccmax2 a*xij^2 + 2*h*xij*yij + b*yij^2]
                
                %                 disp('Valid Result');
                %store valid cut-offs
                rhoc_old = rhoc(i, :);
            else
                %                 disp('Invalid Result');
                vmag(i, :) = NaN;
                vang(i, :) = NaN;
                ar(i, :) = NaN;
                Psi_a(i, :) = NaN;
                Vc(i, :) = NaN;
                args(i, :) = NaN(1, 5);
            end
        end
        i = i + 1;
    end
    
    disp('Finished');
    
    %%
    vctov = Vc ./ vmag;
    vge = vmag .* cos(vang);
    vgn = vmag .* sin(vang);
    vang = vang * 180 / pi;
    vest = [vmag, vang, vge, vgn, ar, Psi_a, vctov];
else
    disp(['maximum correlation value is smaller than hardcoded threshold ', num2str(rho_c)]);
    vest = NaN(1, 7);
    vmag = vest(1);
    vang = vest(2);
    vge = vest(3);
    vgn = vest(4);
    ar = vest(5);
    Psi_a = vest(6);
    vctov = vest(7);
end

%attempt to estimate mean and deviation of drift parameters
try
    ahat = mle(args(~isnan(args(:, 1)), 1));
    hhat = mle(args(~isnan(args(:, 2)), 2));
    bhat = mle(args(~isnan(args(:, 3)), 3));
    fhat = mle(args(~isnan(args(:, 4)), 4));
    ghat = mle(args(~isnan(args(:, 5)), 5));
catch
    ahat = NaN(1, 2);
    hhat = NaN(1, 2);
    bhat = NaN(1, 2);
    fhat = NaN(1, 2);
    ghat = NaN(1, 2);
end
argshat = [ahat, hhat, bhat, fhat, ghat];
try
    vmaghat = mle(vmag(~isnan(vmag)));
catch
    vmaghat = NaN(1, 2);
end
try
    vanghat = mle(vang(~isnan(vang)));
catch
    vanghat = NaN(1, 2);
end
try
    vgehat = mle(vge(~isnan(vge)));
catch
    vgehat = NaN(1, 2);
end
try
    vgnhat = mle(vgn(~isnan(vgn)));
catch
    vgnhat = NaN(1, 2);
end
try
    arhat = mle(ar(~isnan(ar)));
catch
    arhat = NaN(1, 2);
end
try
    psiahat = mle(Psi_a(~isnan(Psi_a)));
catch
    psiahat = NaN(1, 2);
end
try
    vctovhat = mle(vctov(~isnan(vctov)));
catch
    vctovhat = NaN(1, 2);
end

% return;
%attempt to plot drift parameters
set(gcf, 'papersize', [8, 12], ...
    'paperposition', [0, 0, 8, 12], ...
    'paperpositionmode', 'auto', ...
    'position', [0, 0, 8, 12]);
[sp, ~] = tight_subplot(size(vest, 2), 1, [0.02, 0], [0.075, 0.05], [0.15, 0.01]);
for subi = 1:size(vest, 2)
    %     sp(subi) = subplot(size(vest, 2), 1, subi);
    plot(sp(subi), rhoc, vest(:, subi), 'k');
    switch subi
        case 1
            lblstr = ['(a) $v$ [m/s]'];
            %             set(sp(subi),'ytick',0:1000:3000);
            %             ylim([0 3000]);
            %             title(['SAGA estimates as functions of cut-off \rho_c, ',...
            %                 '\Delta{\rho_c} = ',num2str(dt)]);
        case 2
            lblstr = ['(b) $\theta$ [$^\circ$]'];
            set(sp(subi), 'ytick', -180:90:180);
            ylim(sp(subi), [-180, 180]);
        case 3
            lblstr = ['(c) $v_{ge}$'];
        case 4
            lblstr = ['(d) $v_{gn}$'];
        case 5
            lblstr = ['(e) AR'];
            ylim(sp(subi), [0, 10]);
        case 6
            lblstr = ['(f) $\Psi_{\alpha}$ [$^\circ$]'];
            set(sp(subi), 'ytick', -180:90:180);
            ylim(sp(subi), [-180, 180]);
        case 7
            lblstr = ['(g) $v_c / v $'];
    end
    ylabel(sp(subi), lblstr);
    xlim(sp, [rho_c, 1]);
    xtick0 = get(sp(subi), 'xtick');
    %     set(sp(subi),'xtick',sort([min(CCVAL) xtick0 max(CCVAL)]));
    if subi ~= size(vest, 2)
        set(sp(subi), 'XTickLabel', []);
    else
        str2 = ['SAGA estimates vs $\rho_{cutoff}$', ...
            ', $\rho_{min}$ = ', num2str(min(CCVAL), '%0.3f'), ...
            ', $\rho_{{max}}$ = ', num2str(max(CCVAL), '%0.3f')];
        title(sp(1), str2);
        xlabel('$\rho_{cutoff}$');
    end
end
saveas(gcf, '../rhoc.eps', 'epsc2');
close;
end
