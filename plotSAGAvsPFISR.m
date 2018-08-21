function varargout = plotSAGAvsPFISR(prnlist, tstt, vflag)
[prop, op_path0, MEGAVEST_path] = ver_chk;

close all;
warning off;
weimer_path = [MEGAVEST_path, filesep, 'Weimer', filesep, 'runs', filesep];
medium_blue = [0, 0.447, 0.741];
light_blue = [0.301, 0.745, 0.933];

if nargin ~= 0
    dbstop if error;
%     prnlist = [23, 13, 10];
%     tstt = [2013, 12, 8, 4, 3, 0];
    %     prnlist = [23 23];
    %     tstt = [2013, 12, 8, 3, 43, 0];
    %     prnlist = [8, 9, 26];
    %     tstt = [2013, 12, 8, 7, 0, 0];
    %     prnlist = [27, 22, 18];
    %     tstt = [2015, 3, 17, 13, 0, 0];
    %     prnlist = [29];
    %     tstt = [2014, 2, 20, 11, 20, 0];
%     prnlist = [3];
%     tstt = [2015, 10, 7, 6, 2, 0];
    switch length(prnlist)
        case 0
            vflag = 'SAGA';
            args = {[0.03, 0], [0.11, 0.05], [0.11, 0.05]};
            ffflag = 1;
        otherwise
            vflag = 'debug';
            args = {[0.05, 0], [0.18, 0.05], [0.11, 0.05]};
            ffflag = 1:2;
    end
    prnlist = unique(prnlist, 'stable');
    
    year = num2str(tstt(1));
    doy = num2str(floor(datenum(tstt)-datenum([tstt(1), zeros(1, 5)])), '%03i');
    matfilestruct = dir([MEGAVEST_path, 'xcorr_*.mat']);
    
    col_est = [3, 5, 7, 9, 11, 13, 15]; % vmag,vang,vge,vgn,ar,psi_a,vc
    switch vflag
        case 'vmag_vang'
            %plot vgmag and vgang
            col = [3, 5];
            lati = 4;
        case 've_vn'
            col = [7, 9];
            lati = 4;
        case 'ar_psia_vctov'
            col = [11, 13, 15];
            lati = 4;
        case 'SAGA'
            col = [3, 5, 11, 13, 15, 17];
            lati = 4;
        otherwise
            col = [3, 5];
            lati = 4;
    end
    
    for fflag = ffflag
        figure;
        for tau = [60]
            for ilat = lati
                if length(col) > 2
                    set(gcf, 'papersize', [8, 2 * length(col)], ...
                        'paperposition', [0, 0, 8, 2 * length(col)], ...
                        'paperpositionmode', 'auto', ...
                        'position', [0, 0, 8, 2 * length(col)]);
                end
                [sp, ~] = tight_subplot(length(col), 1, args{:});
                for kk = 1:length(prnlist)
                    switch kk
                        case 1
                            lcolor = [1, 0, 1];
                            lcolor2 = 0.66 * lcolor;
                            %                     lcolor2 = [0.5 0.5 0.5]
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
                    end
                    matfilesname = dir([MEGAVEST_path, ...
                        'xcorr_', year, '_', doy, '_PRN', num2str(prnlist(kk)), ...
                        datestr(tstt, '_'), num2str(tau, '*_%is.mat')]);
                    weimermatfile = dir([weimer_path, 'weimer_', year, '_', doy, ...
                        datestr(tstt, '_HHMM'), '*.mat']);
                    try
                        load([weimer_path, weimermatfile.name])
                        %         plus = find(abs(v_weimer(:,3)) == v_weimer(:,3));
                        %         minus = find(abs(v_weimer(:,3)) ~= v_weimer(:,3));
                        %         v_weimer(plus,3) = v_weimer(plus,3) - 180;
                        %         v_weimer(minus,3) = v_weimer(minus,3) + 180;
                    catch
                        %             fprintf('No weimer ExB drifts measurments for this time period yet \n')
                        v_weimer = NaN(1, 5);
                    end
                    if ~isempty(matfilesname)
                        for fi = 1:length(matfilesname)
                            load([MEGAVEST_path, matfilesname(fi).name]);
                            tlimprn(kk, fi, :) = tlim / 24 / 3600 + init_time;
                            
                            switch fflag
                                case 1
                                    fflagstr = 'vc';
                                    filter = imag(ESTV(:, 15)) == 0 & ...
                                        real(ESTV(:, 15)) > 0 & ...
                                        real(ESTV(:, 15)) < ESTV(:, 3);
                                    %                             filter = real(ESTV(:, 15)) > 0 & ...
                                    %                                 real(ESTV(:, 15)) < ESTV(:, 3);
                                case 2
                                    fflagstr = 'ar';
                                    arcap = 20;
                                    percentcap = 80;
                                    filter = ESTV(:, 11) <= arcap & ESTV(:, end) >= percentcap;
                            end
                            
                            %                 filter = ~isnan(ESTV(:, 3));
                            nnan = find(isnan(ESTV(:, 3)));
                            ngood = find(filter);
                            length(ngood)
                            nbad = find(~filter & ~isnan(ESTV(:, 3)));
                            
                            veststats = [];
                            for iiii = 1:length(col_est)
                                veststats(iiii, :) = median(ESTV(~isnan(ESTV(:, col_est(iiii))), col_est(iiii)));
                            end
                            save([MEGAVEST_path, matfilesname(fi).name], 'veststats', '-append');
                            %                 ESTV = ESTV(filter,:);
                            ESTV_ = ESTV';
                            for subi = 1:length(col)
                                if col(subi) == size(ESTV, 2)
                                    vest = ESTV_([1:2, col(subi)], :);
                                    vest = [vest; NaN(1, size(vest, 2))];
                                    vest(end-1, vest(end-1, :) == 0) = 0.999;
                                else
                                    vest = ESTV_([1:2, col(subi), col(subi) + 1], :);
                                    estnotnan{kk, fi} = ngood;
                                end
                                
                                switch col(subi)
                                    case 5
                                        datas{kk, fi} = [mean(vest(1:2, :)); ...
                                            ESTV_([col(subi-1), col(subi-1) + 1], :); ...
                                            vest([end - 1, end], :)]';
                                    case 15
                                        vest(end-1:end, :) = real(vest(end-1:end, :));
                                end
                                tc{kk, fi} = mean(vest(1:2, :))';
                                plotconfig = {marker, 'LineWidth', 1, ...
                                    'Markersize', 4,};
                                %                             hold(sp(subi),'on');
                                h(subi, kk) = errorbar(sp(subi), tc{kk, fi}(ngood), vest(end-1, ngood), vest(end, ngood), ...
                                    plotconfig{:}, 'color', lcolor, 'markerfacecolor', lcolor);
                                %                         sqrt(meansqr(vest(end, ngood) ./ vest(end-1, ngood))) * 100
                                rms((vest(end, ngood) - median(vest(end, ngood)))./vest(end-1, ngood)) * 100
                                hold(sp(subi), 'on');
                                if length(prnlist) == 1
                                    h1(subi, kk) = errorbar(sp(subi), tc{kk, fi}(nbad), vest(end-1, nbad), vest(end, nbad), ...
                                        'x', 'LineWidth', 1, 'Markersize', 6, ...
                                        'color', lcolor2, 'markerfacecolor', lcolor2);
                                end
                                %                     errorbar(sp(subi), tc{kk, fi}(nnan), vest(end-1, nnan), vest(end, nnan), ...
                                %                         plotconfig{:}, 'color', lcolor3, 'markerfacecolor', 0.33*lcolor);
                                %                     plot(TIME,VEST,'LineWidth',2,'color',lcolor);
                                %                     clear TIME VEST;
                                %                     plot_envelope(sp(subi),tc,vest(end-1,:),vest(end,:),lcolor);
                                %                     w = plot(v_weimer(:,1),v_weimer(:,col_w(subi)),'r*-','LineWidth',1);
                                %                             hold(sp(subi),'off');
                                switch col(subi)
                                    case 17
                                        titlestr = ['Estimate Validity Percentage'];
                                        ylabel(sp(subi), '$K/N$ [ \% ]');
                                        %                                     ylim(sp(subi), [0.099, 200]);
                                        set(sp(subi), 'yscale', 'log');
                                        %                                     ylim(sp(subi), [0.099, 200]);
                                    case 3
                                        titlestr = ['Horizontal Drift Magnitude'];
                                        %                                     set(sp(subi), 'ytick', 0:500:2500);
                                        %                                     ylim(sp(subi), [0, 2500]);
                                        ylabel(sp(subi), '$^d v ^{p}$ [m/s]');
                                        set(sp(subi), 'ytick', 0:500:3500);
                                        ylim(sp(subi), [0, 3500]);
                                        %                             set(sp(subi), 'yscale', 'log');
                                    case 5
                                        titlestr = ['Horizontal Drift Orientation Angle'];
                                        set(sp(subi), 'ytick', 0:90:360);
                                        %                             set(sp(subi),'yticklabel',{'West','South','East','North','West'});
                                        ylim(sp(subi), [0, 360]);
                                        ylabel(sp(subi), '$^d \theta ^{p}$ [ $^\circ$ ]');
                                    case 11
                                        titlestr = ['SAGA Axial Ratio of Anisotropy Ellipse'];
                                        set(sp(subi), 'yscale', 'log');
                                        %                                     ylim(sp(subi), [1, 110]);
                                    case 13
                                        titlestr = ['SAGA Anisotropy Ellipse Orientation Angle'];
                                        ylabel(sp(subi), '$\Psi_{\alpha}$ [ $^\circ$ ]');
                                        set(sp(subi), 'ytick', 0:90:180);
                                        ylim(sp(subi), [0, 180]);
                                    case 7
                                        titlestr = ['SAGA Zonal Drift'];
                                        ylim(sp(subi), [-2500, 2500]);
                                        ylabel(sp(subi), '$v_{ge}$ [m/s]');
                                    case 9
                                        titlestr = ['SAGA Meridional Drift'];
                                        ylim(sp(subi), [-2500, 2500]);
                                        ylabel(sp(subi), '$v_{gn}$ [m/s]');
                                    case 15
                                        titlestr = ['SAGA Characteristic Velocity'];
                                        set(sp(subi), 'ytick', 0:500:2500);
                                        ylim(sp(subi), [0, 2500]);
                                        ylabel(sp(subi), '$v_c$ [m/s]');
                                end
                                title(sp(subi), [num2str(subi, '(%i)'), titlestr]);
                                %                             hold(sp(subi),'on');
                            end
                        end
                    else
                        return;
                    end
                end
                tmin = min(tlimprn(:, :, 1));
                tmax = max(tlimprn(:, :, 2));
                if (tmax - tmin) * 24 * 3600 <= 300
                    ticklbl = 'HH:MM:SS';
                    rotang = 0;
                else
                    ticklbl = 'HH:MM';
                    rotang = 15;
                end
                set(sp, 'xlim', [tmin, tmax]);
                for subi = 1:length(col)
                    datetick(sp(subi), 'x', ticklbl);
                    if subi ~= length(col)
                        set(sp(subi), 'XTickLabel', []);
                    else
                        set(sp(subi), 'xticklabelrotation', rotang);
                        xlabel(['Time [HH:MM UT] on: ', datestr(init_time, 'mm/dd/yyyy')]);
                    end
                end
                t_pfisr = datevec([tmin, tmax])
                mflag = 'lat';
                mflag = 'dtau';
                [megadata, lat, dtau, ebar] = plotPFISRvs(t_pfisr(1, :), t_pfisr(2, :), mflag, vflag);
                if ilat ~= 1
                    t_pfisr
                    [lat, ebar(:, :, 1)']
                    %             keyboard;
                end
                if ~isempty(megadata)
                    for jj = 1 %length(dtau)
                        latdata = megadata(:, ilat, jj);
                        for jjj = 1:2
                            tp = latdata{jjj}(:, 1);
                            vp = latdata{jjj}(:, 2);
                            ep = latdata{jjj}(:, 3);
                            for kk = 1:length(prnlist)
                                for fi = 1:length(matfilesname)
                                    vinterp = interp1(tp, vp, tc{kk, fi}, 'linear', 'extrap');
                                    einterp = NaN(size(vinterp));
                                    eb(kk) = errorbar(sp(jjj), tp, vp, ep, ...
                                        '.', 'color', medium_blue, ...
                                        'markerfacecolor', medium_blue, ...
                                        'Markersize', 4);
                                    if isempty(tc{kk, fi}(estnotnan{kk, fi}))
                                        ebi(kk) = errorbar(sp(jjj), tc{kk, fi}, ...
                                            einterp, einterp, ...
                                            's', 'color', light_blue, ...
                                            'markerfacecolor', light_blue, ...
                                            'Markersize', 4);
                                    else
                                        ebi(kk) = errorbar(sp(jjj), tc{kk, fi}(estnotnan{kk, fi}), ...
                                            vinterp(estnotnan{kk, fi}), einterp(estnotnan{kk, fi}), ...
                                            's', 'color', light_blue, ...
                                            'markerfacecolor', light_blue, ...
                                            'Markersize', 4);
                                    end
                                    datai{kk, fi}(:, 1) = tc{kk, fi};
                                    datai{kk, fi}(:, [2 * (jjj - 1) + 2, 2 * (jjj - 1) + 3]) ...
                                        = [vinterp, einterp];
                                end
                            end
                        end
                    end
                    rmsemat = [];
                    
                    for kk = 1:length(prnlist)
                        for fi = 1:length(matfilesname)
                            diffsi = datas{kk, fi} - datai{kk, fi};
                            vdiff = diffsi(estnotnan{kk, fi}, 2);
                            thdiff = diffsi(estnotnan{kk, fi}, 4);
                            tlimprnv = datevec(tlimprn(kk, fi, :));
                            %                     fi, ...
                            %                     tlimprnv(1,4:5), ...
                            %                     tlimprnv(2,4:5), ...
                            rmse = [prnlist(kk), ...
                                length(estnotnan{kk, fi}), ...
                                length(tc{kk, fi}), ...
                                lat(ilat), ...
                                dtau(jj), ...
                                rms(vdiff./datai{kk, fi}(estnotnan{kk, fi}, 2)) * 100, ...
                                rms(thdiff./datai{kk, fi}(estnotnan{kk, fi}, 4)) * 100];
                            rmsemat = [rmsemat; rmse];
                        end
                    end
                    if length(prnlist) == 1
                        lgs = [h(end, :), h1(end, :), eb(kk), ebi(kk)];
                    else
                        lgs = [h(end, :), eb(kk), ebi(kk)];
                    end
                    latstr = num2str(lat(ilat));
                else
                    lgs = [h(end, :), h1(end, :)];
                    latstr = '';
                end
                if length(prnlist) == 1
                    liststr = {num2str(prnlist', 'PRN%i$\\surd$'); ...
                        num2str(prnlist', 'PRN%i$\\times$')};
                else
                    liststr = cellstr(num2str(prnlist', 'PRN%i'));
                end
                lg = legend(lgs, ...
                    [liststr; ...
                    ['PFISR ', latstr, '$^\circ$']; ...
                    ['PFISR ', latstr, '$^\circ$Interp']], ...
                    'location', 'North', 'orientation', 'horizontal');
                lgpos = get(lg, 'Position');
                spposu = get(sp(end-1), 'outerPosition');
                spposl = get(sp(end), 'outerPosition');
                height = spposu(2) - (spposl(2) + spposl(4));
                set(lg, 'position', ...
                    [(1 - lgpos(3)) / 2, 0.005, lgpos(3), lgpos(4)]);
                tmp = strjoin(cellstr(num2str(prnlist', '%02i'))', '_');
                a = {'SAGA', vflag, 'PRN', tmp, year, doy, datestr(tstt, 'HHMMUT'), ...
                    num2str(0.65), num2str(tau, '%gs'), latstr, fflagstr};
                plotfilename = strjoin(a, '_');
%                 saveas(gcf, [op_path0, plotfilename, '.eps'], 'epsc2');
                saveas(gcf, [op_path0, plotfilename, '.png'], 'png');
                %         try
                %             print([op_path0, plotfilename,'.pdf'], '-dpdf', '-fillpage');
                %         catch
                %             print([op_path0, plotfilename,'.pdf'], '-dpdf');
                %         end
                close;
                clear tc datai;
                if strcmp(vflag, 'debug')
                    varargout = {rmsemat};
                    disp(rmsemat);
                else
                    varargout = {};
                end
            end
        end
    end
end
