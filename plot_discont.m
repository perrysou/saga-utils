function [tspan] = plot_discont(times, dt)
% times = [1 3 4 6 7 10 11 13]';
% dt = 2;
% discont_proc(times,times,dt)
tspan = [];
j = 1;
while j <= length(times)
    idx = find(times-times(j) <= dt & times-times(j) >= 0);
    %             if ~isempty(idx)
    %                 tt = [times(j);times(idx)];
    %                 tt = [tt(1);tt(end)];
    %                 j = idx(end)+1;
    %             else
    %                 tt = [times(j);times(j)];
    %                 j = j+1;
    %             end
    tt = [times(j); times(idx)];
    tt = [tt(1); tt(end)];
    j = idx(end) + 1;
    %     tspan = [tspan;tt;NaN;NaN];
    tspan = [tspan; tt];
end
odd = tspan(1:2:end);
even = tspan(2:2:end);
tspan_new = [odd, even];
% plot(tspan,'-');
% datetick('y','HH:MM','keeplimits');
close;
end
