function [] = plot_el(DATA_el, kk, dt, color, tlim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sv_el = find(DATA_el(:, 4) == kk);
if ~isempty(sv_el)
    ortw_el = DATA_el(sv_el, 1);
    orts_el = DATA_el(sv_el, 2);
    el = DATA_el(sv_el, 3);
    time_el = gps2utc([ortw_el, orts_el], 0);
    time_el = datenum(time_el);
    [time_el_e, el_e] = discont_proc(time_el, el, dt);
    plot(time_el_e, el_e/90, 'LineWidth', 1, 'Color', color);
    set(gca, 'FontSize', 8, 'YTick', 1/90*[30, 90], ... ,
        'YTickLabel', {'30,1/3', '90,1'}, 'YGrid', 'on');
    ylim([0, 1.1]);
    xlabel(['PRN', num2str(kk)]);
end

if isempty(tlim)
    datetick('x', 'HH', 'keeplimits');
else
    xlim(tlim);
    datetick('x', 'HH', 'keeplimits');
end
hold on;

end
