function [] = plotlz(prnlist, tstt, vflag)
[prop, op_path0, MEGAVEST_path] = ver_chk;

close all;
dbstop if error;
warning off;
weimer_path = [MEGAVEST_path, filesep, 'Weimer', filesep, 'runs', filesep];
medium_blue = [0, 0.447, 0.741];
light_blue = [0.301, 0.745, 0.933];
format = 'png';
% format = 'epsc2';

if nargin == 0
    dbstop if error;
    %     vflag = 'vmag_vang';
    %     vflag = 've_vn';
    %     vflag = 'ar_psia_vctov';
    vflag = 'debug';
    prnlist = [23, 13, 10];
    tstt = [2013, 12, 8, 4, 3, 0];
    prnlist = [23];
    tstt = [2013, 12, 8, 3, 43, 0];
    %     prnlist = [8, 9, 26];
    %     tstt = [2013, 12, 8, 7, 0, 0];
    %     prnlist = [25, 29, 31];
    %     tstt = [2013, 12, 8, 16, 0, 0];
%         prnlist = [29];
%         tstt = [2014, 2, 20, 11, 20, 0];
end
year = num2str(tstt(1));
doy = num2str(floor(datenum(tstt)-datenum([tstt(1), zeros(1, 5)])), '%03i');
matfilestruct = dir([MEGAVEST_path, 'xcorr_*.mat']);
tau = 60;
col = [3 4 5];
figbox = figure;
[sp, ~] = tight_subplot(length(col), 1, [0.05, 0], [0.18, 0.05], [0.11, 0.05]);
matfilesname = dir([MEGAVEST_path, ...
    'xcorr_', year, '_', doy, '_PRN', num2str(prnlist), ...
    datestr(tstt, '_'), num2str(tau,'*_%is.mat')]);

    if ~isempty(matfilesname)
        for fi = 1:length(matfilesname)

            load([MEGAVEST_path, matfilesname(fi).name]);            
            for kk = 1:length(sitenum_op)
                switch kk
                    case 1
                        lcolor = [1, 0, 1];
                        lcolor2 = 0.66 * lcolor;
                        lcolor3 = 0.33 * lcolor;
                        marker = 'd';
                    case 2
                        lcolor = [1, 0, 0];
                        lcolor2 = 0.66 * lcolor;
                        lcolor3 = 0.33 * lcolor;
                        marker = '^';
                    case 3
                        lcolor = [0, 0, 0];
                        lcolor2 = 0.5 + lcolor;
                        lcolor3 = 0.75 + lcolor;
                        marker = 'o';
                    case 4
                        lcolor = [0, 1, 0];
                        lcolor2 = 0.66 * lcolor;
                        lcolor3 = 0.33 * lcolor;
                        marker = 's';
                    case 5
                        lcolor = [0, 0, 1];
                        lcolor2 = 0.66 * lcolor;
                        lcolor3 = 0.33 * lcolor;
                        marker = '>';
                end
        
                tlimprn(fi,:) = tlim / 24 / 3600 + init_time;

                
                filter = MEGA_LZ(:, kk, 5) >= 0;
                ngood = find(filter);
                nbad = find(~filter & ~isnan(MEGA_LZ(:, kk, 3)));

                for subi = 1:length(col)
                    vest = MEGA_LZ(:,kk,[1:2, col(subi)]);
                    tc{kk, fi} = mean(vest(:,:,1:2),3);
                    t0{kk, fi} = vest(:,:,1);
                    tf{kk, fi} = vest(:,:,2);
                    plotconfig = {marker, 'LineWidth', 1, ...
                        'Markersize', 4,};
%                     h(subi, kk) = errorbar(sp(subi), tc{kk, fi}(ngood), vest(ngood,:,end), [], ...
%                         plotconfig{:}, 'color', lcolor, 'markerfacecolor', lcolor);
%                     hold(sp(subi), 'on');

                    
                    x_kk = MEGA_LZ(:, kk, 1);
                    y_kk = MEGA_LZ(:, kk, 2);
                    bottomheight_kk = MEGA_LZ(:, kk, 4) - MEGA_LZ(:, kk, 3); 
                    topheight_kk = MEGA_LZ(:, kk, 4);
                    
                    if col(subi) == 5
                        rect(kk, :) = patch(sp(subi), ...
                            [x_kk, y_kk, y_kk, x_kk]', ...
                            [bottomheight_kk, bottomheight_kk, topheight_kk, topheight_kk]', ...
                            'k', 'edgecolor', rx_color(sitenum_op{kk}), ...
                            'linewidth', 1.5, 'facealpha', 1 / length(sitenum_op));
                        set(sp(subi), 'layer', 'top');
                    end   
                    
                    x = MEGA_LZ(:, :, 1);
                    y = MEGA_LZ(:, :, 2);
                    if col(subi) == 3
%                         bottomheight = MEGA_LZ(:, :, 4) - MEGA_LZ(:, :, 3); 
                        thickness = MEGA_LZ(:, :, col(subi));  
                        plotBox(sp(subi), x(:, 1), y(:, 1), thickness);
                    end
                    
                    if col(subi) == 4
                        % DIY version of boxplot
                        topheight = MEGA_LZ(:, :, col(subi));                        
                        plotBox(sp(subi), x(:, 1), y(:, 1), topheight);
                    end                    
                        
                        % debugging with the matlab boxplot
%                         boxplot(sp(subi + 1), bottomheight');
%                     if length(prnlist) == 1
%                         h1(subi, kk) = errorbar(sp(subi), tc{kk, fi}(nbad), vest(nbad,:,end), [], ...
%                             plotconfig{:}, 'color', lcolor2, 'markerfacecolor', lcolor2);
%                     end
                    switch col(subi)
                        case 3
                            titlestr = ['SAGA Thickness'];
                            ylabel(sp(subi), '$\hat{L}$ [km]');
                        case 4
                            titlestr = ['SAGA Top Height'];
                            ylabel(sp(subi), '$\hat{z}$ [km]');
                        case 5
%                             titlestr = ['SAGA Optimal Fit Error'];
%                             ylabel(sp(subi), '$\epsilon^2_{min}$');
                            titlestr = ['SAGA Estimated Irregularity Layer'];
                            ylabel(sp(subi), 'Altitude [km]');
                    end
                    title(sp(subi), [num2str(subi, '(%i)'), titlestr]);
                    %                             hold(sp(subi),'on');
                end
            end
        end
        set(sp, 'ytick', 100:100:650, 'yticklabel', 100:100:650);
        ylim(sp, [50, 650]);
    else
        return;
    end

tmin = min(tlimprn(:, 1));
tmax = max(tlimprn(:, 2));
if (tmax - tmin) * 24 * 3600 <= 300
    ticklbl = 'HH:MM:SS';
    rotang = 0;
else
    ticklbl = 'HH:MM';
    rotang = 15;
end
set(sp, 'xlim', [tmin, tmax]);
for subi = 1:length(col)
    datetick(sp(subi), 'x', ticklbl, 'keeplimits');
    if subi ~= length(col)
        set(sp(subi), 'XTickLabel', []);
    else
%         set(sp(subi), 'xticklabelrotation', rotang);
        xlabel(['Time [HH:MM UT] on: ', datestr(init_time, 'mm/dd/yyyy')]);
    end
end

% adjust the legend
lg = legend(rect, sitenum_op, ...
    'location','North','orientation','horizontal');
lgpos = get(lg, 'Position');
spposu = get(sp(end-1), 'outerPosition');
spposl = get(sp(end), 'outerPosition');
height = spposu(2) - (spposl(2) + spposl(4));
set(lg, 'position', ...
    [(1 - lgpos(3)) / 2, 0.005, lgpos(3), lgpos(4)]);

h = get(sp(end), 'children');
figcomp = figure;
[~, ~, ~, op_path] = plotPFISR_NeTe(tmin, tmax, 'Ne');
copyobj(h, gca);
set(gca, 'xtick', get(sp(end), 'xtick'));
datetick('x', 'HH:MM', 'keepticks');
set(gca, 'xlim', [tmin, tmax]);

saveas(figbox,[op_path,'SAGA_LZ_PRN',num2str(prnlist),'_',year,'_',doy], format);
saveas(figcomp, [op_path,'PFISR_SAGA_LZ_PRN',num2str(prnlist),'_',year,'_',doy], format);

close all;
end