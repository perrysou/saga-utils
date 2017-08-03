function plotellipse(AR, ma, mi, Psi_a, V, Psi_v, cut_off, varargin)
% plot anisotropy ellipse of the ground pattern

if nargin == 0
    load('ellipse.mat');
else
    save('ellipse.mat')
end

[CCVAL, COMBOS, RHO0, BASELINE, RCVRID] = varargin{:};
SITEID = rx2site(RCVRID);
XII = BASELINE(:, 1); YII = BASELINE(:, 2);
XJJ = BASELINE(:, 3); YJJ = BASELINE(:, 4);
% mean(a) .* XIJ.^2 + 2.* mean(h) .* XIJ .* YIJ + mean(b) .* YIJ.^2 - max(CCVAL,[],2)
uniqueij = unique(COMBOS, 'rows', 'stable');
cmap = colormap(hsv(length(uniqueij)));
uniquerho0 = unique(RHO0, 'rows', 'stable');
if ~isnan(V)
    for j = 1:size(uniqueij, 1)
        ijrows = find(COMBOS(:, 1) == uniqueij(j, 1) & COMBOS(:, 2) == uniqueij(j, 2));
        offset = [unique(XII(ijrows)); unique(YII(ijrows))];
        edge = [unique(XJJ(ijrows)); unique(YJJ(ijrows))];
        for ccval = uniquerho0(j,2)
            %ellipse resolution N
            N = 100;
            th = linspace(0, 2*pi, N);
            % unrotated xaxis and yaxis
            x = mi * AR * cos(th) * sqrt(1/ccval) / 1;
            y = mi * sin(th) * sqrt(1/ccval) / 1;
            R = [cos(Psi_a), -sin(Psi_a); ...
                sin(Psi_a), cos(Psi_a)];
            % rotated coords
            xyrot = R * [x; y];
            xrot = offset(1) + xyrot(1,:);
            yrot = offset(2) + xyrot(2,:);
            
            major = R * [max(x), min(x); 0, 0];
            minor = R * [0, 0; max(y), min(y)];
            %plot ellipse
            hold on;
            hellipse(j) = plot(xrot, yrot, 'color', cmap(j,:), 'linewidth', 3);
            %                 drawnow;
            plot(offset(1)+major(1,:), ...
                offset(2)+major(2,:), ':k');
            plot(offset(1)+minor(1,:), ...
                offset(2)+minor(2,:), ':k');
        end
        %         plot([unique(XII(ijrows)) unique(XJJ(ijrows))], ...
        %             [unique(YII(ijrows)) unique(YJJ(ijrows))],'o-');
        plot(offset(1), offset(2), 'o', 'Color', ...
            rx_color(RCVRID(uniqueij(j, 1),:)), ...
            'MarkerFaceColor', rx_color((RCVRID(uniqueij(j, 1),:))));
        plot(edge(1), edge(2), 'o', 'Color', ...
            rx_color(RCVRID(uniqueij(j, 2),:)), ...
            'MarkerFaceColor', rx_color((RCVRID(uniqueij(j, 2),:))));
        text(offset(1), offset(2), SITEID{uniqueij(j, 1)});
        text(edge(1), edge(2), SITEID{uniqueij(j, 2)});
        lg{j} = [SITEID{uniqueij(j, 1)}, '\&', SITEID{uniqueij(j, 2)}, ...
            ',$\rho(0)=$', num2str(ccval)];
        if j == 1
            plot([offset(1), cos(Psi_v) * V], ...
                [offset(2), sin(Psi_v) * V], 'k', 'linewidth', 2);
        end
    end
    legend(hellipse, lg);
end

grid on; box on; axis equal;
title({'Anisotropy Ellipse of the Ground Diffraction Pattern'; ...
    ['$\rho_{c} = ', num2str(cut_off), ', \|V\| = ', ...
    num2str(V, '%.2f'), 'm/s, \theta = ', num2str(Psi_v*180/pi, '%.2f'), ...
    '^\circ', ', AR = ', num2str(AR, '%.2f'), ', \Psi_a = ', ...
    num2str(Psi_a*180/pi, '%.2f'), '^\circ$'];});
end


