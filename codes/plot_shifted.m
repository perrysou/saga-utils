function [ph_truncated, pwr_truncated] = ...
    plot_shifted(xdata, tpeak, rcvr_op, sitenum_op, combos,  Fs)
% load('lz.mat');
combos_auto = [1:size(rcvr_op, 1); 1:size(rcvr_op, 1)]';
[~, index] = sortrows([combos; fliplr(combos); combos_auto]);
tpeaks = [tpeak'; -tpeak'; zeros(size(rcvr_op, 1), 1)];
tpeaks_new = tpeaks(index);
tpeaks_array = reshape(tpeaks_new, [], size(rcvr_op, 1));

for rr = 1:size(rcvr_op, 1)
    tpeaks_rr = tpeaks_array(:, rr);
    leadnum(rr) = length(find(tpeaks_rr < 0));
    lagnum(rr) = size(rcvr_op, 1) - 1 - leadnum(rr);
end
rr_lead = find(leadnum == max(leadnum));
rr_lag = find(lagnum == max(lagnum));

% tmat(:,1) = tmat(:,1) - 2.79;
% tmat(:,2) = tmat(:,2) - 2.53;
% tmat(:,3) = tmat(:,3) - 1.62;
% tmat(:,4) = tmat(:,4);
% tmat(:,5) = tmat(:,5) - 2.13;

fig_all = figure;
[ph_shifted, pwr_shifted, ph_truncated, pwr_truncated] = deal(cell(1, size(rcvr_op, 1)));
for rr = 1:size(rcvr_op, 1)
    ph_shifted{rr} = xdata{rr}(tpeaks_array(rr_lead, rr) * Fs + 1 : end, 3);
    pwr_shifted{rr} = xdata{rr}(tpeaks_array(rr_lead, rr) * Fs + 1 : end, 2);
end

[sp, ~] = tight_subplot(2, 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
for rr = 1:size(rcvr_op, 1)
    ph_truncated{rr} = ph_shifted{rr}(1 : length(ph_shifted{rr_lag}));  
    pwr_truncated{rr} = pwr_shifted{rr}(1 : length(pwr_shifted{rr_lag})); 
    plot(sp(2), ph_shifted{rr}, 'color', rx_color(rcvr_op(rr,:)));
    plot(sp(1), pwr_shifted{rr}, 'color', rx_color(rcvr_op(rr,:)));
    hold(sp(1), 'on'); hold(sp(end), 'on');
end
set(sp(1 : (end - 1)), 'xticklabel', []);
legend(sp(1), sitenum_op, 'orientation', 'horizontal');
tightfig;
saveas(fig_all, '../../../all', 'png');
close;

fig_each = figure;
if size(rcvr_op, 1) > 2
    set(gcf, 'papersize', [8, 2*size(rcvr_op, 1)], ...
        'paperposition', [0, 0, 8, 2*size(rcvr_op, 1)], ...
        'paperpositionmode', 'auto', ...
        'position', [0, 0, 8, 2*size(rcvr_op, 1)]);
end
[sp, ~] = tight_subplot(size(rcvr_op, 1), 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
for rr = 1:size(rcvr_op, 1)  
%     plot(sp(rr), ph_shifted{rr});
    plot(sp(rr), ph_truncated{rr}, 'color', rx_color(rcvr_op(rr,:)));
    legend(sp(rr), sitenum_op{rr});
end
set(sp(1 : end - 1), 'xticklabel', []);
set(sp(2: 2: end), 'yticklabel', []);
tightfig;
saveas(fig_each, '../../../each', 'png');
close;
