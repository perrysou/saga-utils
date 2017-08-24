function [] = SAGAvsPFISR_plotter(prnlist,tstt,vflag)
clear all; close all;
if verLessThan('matlab','8.4.0') && strcmp('GLNXA64',computer)
    prop = 'clim';
    op_path0 = '/data1/home/ysu27/Dropbox/';
%     op_path0 = '../'
%     op_path0 = '/data1/home/ysu27/Dropbox/Comprehensive_Exam/figs/';
%     op_path0 = 'D:\Dropbox\research-misc\[MTSSP]\2015_Yang_MTSSP\v6\figs\';
    MEGAVEST_path = '/data1/home/ysu27/Dropbox/research/';
elseif verLessThan('matlab','8.4.0') && strcmp('PCWIN64',computer)
    prop = 'clim';
    op_path0 = 'C:\Users\Yang Su\Dropbox\';
%     op_path0 = '../'
%     op_path0 = '/data1/home/ysu27/Dropbox/Comprehensive_Exam/figs/';
%     op_path0 = 'D:\Dropbox\research-misc\[MTSSP]\2015_Yang_MTSSP\v6\figs\';
    MEGAVEST_path = 'C:\Users\Yang Su\Dropbox\research\';
else
    prop = 'limits';
    op_path0 = '/Users/yangsu/Dropbox/';
    op_path0 = '../';
%     op_path0 = '/Users/yangsu/Dropbox/Comprehensive_Exam/figs/';
%     op_path0 = 'D:\Dropbox\research-misc\[MTSSP]\2015_Yang_MTSSP\v6\figs\';
    MEGAVEST_path = '/Users/yangsu/Dropbox/research/';
%     MEGAVEST_path = 'D:\Dropbox\research\';
    MEGAVEST_path = './';
end

figure;
set(gcf,'paperpositionmode','auto');
set(gcf,'units','normalized');

vflag = 'vmag_vang';
% prnlist = [23];
% tstt = [2013 12 8 3 42 0];
% year = num2str(tstt(1));
% doy = num2str(floor(datenum(tstt) - datenum([tstt(1) zeros(1,5)])),'%03i')
% prnlist = [29];
% year = '2014'
% doy = '051'
% prnlist = [27 22];
% year  = '2015'
% doy = '076';
matfilestruct = dir([MEGAVEST_path,'xcorr_*.mat'])

col_est = [3 5 7 9 11 13 15]; % vmag,vang,vge,vgn,ar,psi_a,vc
switch vflag
    case 'vmag_vang'
        %plot vgmag and vgang
        col = [3 5];
    case 've_vn'
        %plot vge and vgn [14 15]
        col = [7 9];
end

for ccmin = 0.6
for ii = 1:50:size(matfilestruct,1)
% for kk = 1:length(prnlist)
    switch mod(ii,6)
        case 0
            lcolor = 'y';
            lcolor = 'k';
            marker = 'o';
        case 1
            lcolor = 'k';
            marker = 'd';
        case 2
            lcolor = 'r';
            lcolor = 'k';
            marker = '^';
        case 3
            lcolor = 'm';
            lcolor = 'k';
            marker = 's';
        case 4
            lcolor = 'b';
            lcolor = 'k';
            marker = 'x';
        case 5
            lcolor = 'g';
            lcolor = 'k';
            marker = 'v';
    end
    matfilesname = matfilestruct(ii)
    if ~isempty(matfilesname)
        load([MEGAVEST_path,matfilesname.name]);        
        for iiii = 1:length(col_est)
            veststats(iiii,:) = median(ESTV(~isnan(ESTV(:,col_est(iiii))),col_est(iiii)));      
        end
        save([MEGAVEST_path,matfilesname.name],'veststats','-append');
       
        for i = find(v_ccmin==ccmin)
            ccmin_str = ['\rho_c = ',num2str(v_ccmin(i))];
            for j = find(v_dtau==60) %dtau  = 60
                dtau_str = ['\Delta\tau = ',num2str(v_dtau(j)),' s'];
                for subi = 1:length(col)
                    sp(subi) = subplot(length(col),1,subi); 
                    box on;
                    switch col(subi)
                        case 3
                            titlestr = ['(a) SAGA Horizontal Drift Magnitude'];
                            set(sp(subi),'ytick',0:500:2500);
                            ylim([0 2500]);
                            ylabel('v_{mag} [m/s]');
                        case 5
                            titlestr = ['(b) SAGA Horizontal Drift Orientation'];
                            set(sp(subi),'ytick',-360:90:360);
%                             set(sp(subi),'yticklabel',{'West','South','East','North','West'});
                            ylim([-360 360]);
                            ylabel('v_{ang} [deg]');
                        case 11
                            titlestr = ['(d) SAGA Axial Ratio of Anisotropy Ellipse'];
                        case 13
                            titlestr = ['(e) SAGA Anisotropy Ellipse Orientation'];
                            ylabel('\Psi_{\alpha} [\circ]');
                            set(sp(subi),'ytick',0:90:180);
                            ylim([0 180]);
                        case 7
                            titlestr = ['(a) SAGA Zonal Drift'];
                            ylim([-2500 2500]);
                            ylabel('v_{ge} [m/s]');
                        case 9
                            titlestr = ['(b) SAGA Meridional Drift'];
                            ylim([-2500 2500]);
                            ylabel('v_{gn} [m/s]');
                        case 15
                            titlestr = ['(c) SAGA Characteristic Velocity v_{char}/v_{mag}'];
                            set(sp(subi),'ytick',0:0.25:1);
                            ylim([0 1]);
                    end  
%                     xlabel([lblstr,ccmin_str,', ',dtau_str]);
                    title(titlestr);
                    ESTV_ = ESTV';
                                       
                    vest = ESTV_([1:2 col(subi) col(subi)+1],:);
                    if col(subi) == 13
                        minus_ind = find(vest(end-1,:)<0);
                        vest(end-1,minus_ind) = 180 + vest(end-1,minus_ind);
                    end
                    tc = mean(vest(1:2,:));
                    init_time_v = datevec(init_time);
                    doy = floor(datenum(init_time_v) - datenum([init_time_v(1) zeros(1,5)]));
                    tc_v = datevec(tc);
                    tc = datenum([repmat([2013 12 8],size(tc_v,1),1) tc_v(:,4:6)]);
%                     time = vest(1:2,:);
%                     for iii = 1:size(vest,1)
%                         TIME(3*(iii-1)+1:3*iii,:) = [time(1,iii);time(2,iii);NaN]; 
%                         VEST(3*(iii-1)+1:3*iii,:) = [vest(3,iii);vest(3,iii);NaN];
%                     end
                    hold on;grid on;
%                     vest(:,17:20)
                    h(subi,ii) = errorbar(sp(subi),tc,vest(end-1,:),vest(end,:),...
                        marker,'color',lcolor,'Markersize',1,'markerfacecolor',lcolor);
                    plot_envelope(sp(subi),tc,vest(end-1,:)',vest(end,:)');
%                     plot(TIME,VEST,'LineWidth',2,'color',lcolor);
%                     clear TIME VEST;
%                     datetick(sp(subi),'x');                         

%                     xlim(tlim/24/3600+init_time);
                end
                
                t_pfisr = datevec(tlim/24/3600+init_time);
                PFISR_subplotter(sp,t_pfisr(1,:),t_pfisr(2,:),'dtau',vflag);
            end
        end
    else
        return;
    end
end

for subi = 1:length(col)
    datetick(sp(subi),'x','HH');
    if subi ~= length(col)
        set(sp(subi),'XTickLabel',[]);
    else
%         xlabel(['Time [HH:MM UT] on: ',datestr(init_time,'mm/dd/yyyy')]);
        xlabel(['Time [HH:MM UT]']);
    end
end
keyboard;


for iii = length(col):-1:1
    op{iii} = get(sp(iii),'outerposition');
    op{iii};
    p{iii} = get(sp(iii),'position');
    p{iii};
end
% %left bottom width height
% set(sp(5),'outerposition',[0    0    1.0000    0.24])
% set(sp(4),'outerposition',[0    0.18    1.0000    0.2])
% set(sp(3),'outerposition',[0    0.33    1.0000    0.2])
% set(sp(2),'outerposition',[0    0.48    1.0000    0.23])
% set(sp(1),'outerposition',[0    0.66   1.0000    0.25])

% lg = legend([h(1,:) eb],[cellstr([repmat('SAGA PRN',length(prnlist),1),num2str(prnlist')]);['PFISR ',lgstr]],...
%     'location','North','orientation','horizontal');
% lg = legend([h(end,:) eb],[cellstr([repmat('SAGA PRN',length(prnlist),1),num2str(prnlist')]);['PFISR ',lgstr]],...
%     'location','South','orientation','horizontal');
%     lblstr = 'Time UT HH:MM')
% get(lg,'Position');
% set(lg,'Position',[0.25 0.95 0.5 0.1]);
% set(lg,'Position',[0.285 0.52 0.4591 0.0356]);

set(findall(gcf,'-property','FontSize'),'FontSize',14);
% set(findall(gcf,'-property','FontSize'),'FontName','Times');
saveas(gcf,[op_path0,'SAGA_',vflag,'_PRN',num2str(prnlist),'_',year,'_',doy,'_',datestr(tstt,'HHMMUT_'),num2str(v_ccmin(i)),'_',num2str(v_dtau(j)),'s'...
'.eps'],'epsc2');

% print('-depsc2',[op_path0,'SAGA_PRNvs_',num2str(v_ccmin(i)),'_',num2str(v_dtau(j)),'s'...
%     '.eps'])
        close;

end
end

function [] = PFISR_subplotter(sp,t1,t2,mflag,vflag)
% t_pfisr = datevec(get(sp(subi),'xlim'));
% mflag = 'lat';
% mflag = 'dtau';
[megadata,lat,dtau] = plotPFISRvs(t1,t2,mflag,vflag);
if ~isempty(megadata)
switch mflag
    case 'dtau'
        for jj = 1
%         for jj = 1:length(dtau)
            latdata = megadata(:,5,jj);
            for jjj = 1:2
                ttt = latdata{jjj}(:,1);
                ttt_v = datevec(ttt);
                ttt = datenum([repmat([2013 12 8],size(ttt_v,1),1) ttt_v(:,4:6)]);
                eb = errorbar(sp(jjj),ttt,latdata{jjj}(:,2),latdata{jjj}(:,3),'.','color','c','Markersize',9,'LineWidth',1);
            end
            drawnow;
            lgstr = [num2str(lat(5)),'\circ'] ;
        end
    case 'lat'
        for jj = 1:length(lat)
            latdata = megadata(:,jj,1);
            for jjj = 1:2                
                ttt = latdata{jjj}(:,1);
                ttt_v = datevec(ttt);
                ttt = datenum([repmat([2013 12 8],size(ttt_v,1),1) ttt_v(:,4:6)]);
                eb = errorbar(sp(jjj),ttt,latdata{jjj}(:,2),latdata{jjj}(:,3),'.','color','c','Markersize',9,'LineWidth',1);
            end
            drawnow;
            lgstr = [num2str(dtau(1)),' s'] ;
        end
end
end
end

