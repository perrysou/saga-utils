function [events] = findhrtimes(year, doy)

dbstop if error;

[prop, op_path, mat_path] = ver_chk;

%number of simultaneously scintillating satellites
% nscint = 3;

if nargin == 0
    year = '2013';
    doy = '342';
%     year = '2014';
%     doy = '051';
%     year = '2015';
%     doy = '076';
end
load([mat_path, 'lrtimes_', year, '_', doy, '.mat']);

% prn = TSP_hr(:,1);
% t0 = TSP_hr(:,2);
% tf = TSP_hr(:,3);
% meansp = TSP_hr(:,4);
% npoints = TSP_hr(:,5);

tlim1 = datenum([str2double(year), 0, 0, 0, 0, 0]) + str2double(doy);
tlim = [tlim1, tlim1 + 1];

clim = [0, 2 * pi];
datalim = clim(2);
cblabel = '$\sigma_\Phi$ [rad]';

prnlist = unique(TSP_hr0(:, 1));


tseg = cell(1, length(prnlist));
for kk = 1:size(prnlist, 1)
    prn = prnlist(kk,:);
    t0 = TSP_hr0(TSP_hr0(:, 1) == prn, 2);
    tf = TSP_hr0(TSP_hr0(:, 1) == prn, 3);
    weight = TSP_hr0(TSP_hr0(:, 1) == prn, 5);
    
    colormap(jet);
    caxis(clim);
    for i = 1:length(t0)
        ind = find(MSP(:, 1) <= tf(i) & MSP(:, 1) >= t0(i) & MSP(:, 3) == prn);
        %         (tf(i) - t0(i)) * 24 * 60;
        h_scatter = scatter(MSP(ind, 1), MSP(ind, 2)/datalim+kk-1, [], MSP(ind, 2), '.');
        hold on;
        %         h_scatter = stairs(MSP(ind,1),MSP(ind,2)/datalim+kk);
    end
    
    sorted = sortrows([t0, tf, weight], 1);
    t0 = sorted(:, 1)';
    tf = sorted(:, 2)';
    weight = sorted(:, 3) / datalim + kk - 1;
    
    tseg{kk} = reshape([t0; tf], [], 1)';
    
    time = reshape([t0; tf; NaN(1, length(t0))], [], 1);
    h_mean = plot(time, reshape(repmat(weight', [3, 1]), [], 1), 'k', 'linewidth', 1);
    %     text(t0,weight,num2str(prn),'VerticalAlignment','Baseline','HorizontalAlignment',...
    %             'Left');
%     text(t0(1), kk-1, num2str(prn), 'VerticalAlignment', 'Baseline', ...
%         'HorizontalAlignment', 'right');
    
end

[T, C] = find_common_times_v2(tseg);
pairs = 1:size(T,1);
events = [];
for row = pairs
    validind(row) = ~isempty(vertcat(T{row,:})); 
    npair(row) = size(C{row,:},2);
end

maxvalid = min(pairs(validind));
for iii = pairs(validind)   
    for jjj = 1:size(T, 2)
        if ~isempty(T{iii, jjj}) 
            t = T{iii, jjj};
            prn_scint = prnlist(C{iii,:}(jjj,:));
            nscint = length(prn_scint);
            tcommon = reshape(t, 2, [])';
            t00 = tcommon(:, 1);
            tff = tcommon(:, 2);
            for i = 1:length(t00)
                ind = find(MSP(:, 1) <= tff(i) & MSP(:, 1) >= t00(i));
                h_scatter = scatter( ...
                    MSP(ind, 1), MSP(ind, 2)/datalim-1+maxvalid-iii, [], MSP(ind, 2), '.');
                hold on;
%                 if iii == maxvalid  
%                     text(tff(i), -1+maxvalid-iii, num2str(prn_scint','%02i,'), ...
%                         'VerticalAlignment', 'Baseline', ...
%                         'HorizontalAlignment', 'left', ...
%                         'edgecolor',[0 0 0]);
%                 end
                if iii == maxvalid 
                    rectangle('position',[t00(i),...
                        find(prnlist==min(prn_scint))-1+maxvalid-iii,...
                        tff(i)-t00(i),...
                        find(prnlist==max(prn_scint))-find(prnlist==min(prn_scint))+1]);
                end
            end
            t_overlap = reshape([reshape(t, 2, []); NaN(1, size(t, 1)/2)], [], 1);
            events = [events; ...
                prn_scint, repmat(t_overlap(1:2)', nscint, 1), ...
                repmat(nscint, nscint, 1)];
            plot(t_overlap, ...
                repmat(mean(MSP(ind, 2))/datalim-1+maxvalid-iii, size(t_overlap)), 'k');
        end
    end
end

events = sortrows(events, [-4, 2]);
events = events(events(:, end) == max(events(:, 4)),:);

legend([h_scatter, h_mean], {'value', 'mean'});
cb = colorbar;
set(cb, prop, clim);
label = get(cb, 'Label');
set(label, 'String', cblabel,'interpreter','latex');

xlim(tlim);
datetick('x', 'HH', 'keeplimits');
pairlabel = num2str(flip(npair(validind)'),'$N_{set}$ = %i');
prnlabel = num2str(prnlist,'PRN%02i');
a = {pairlabel;prnlabel};
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'Ytick', (maxvalid-kk):kk-1, ...
    'YTickLabel', a);
Tick = get(gca, 'xtick');
xlabel(['Time [HH UT] on ', datestr(Tick(1), 'mm/dd/yyyy')]);
title(['Potential Scintillating Times after $\sigma_{\Phi}$ threshold [rad]: ', num2str(spth_hr)]);
plotname = ['hrtimes_', year, '_', doy];
plotpath = [op_path, plotname, '.eps'];
saveas(gcf, plotpath, 'epsc2');
close;


