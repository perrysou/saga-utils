function plot_vel_v3(flag,dtau,ccmin,year,doy,prnlist)
%
if verLessThan('matlab','8.4.0')
    prop = 'clim';
    op_path = '/data1/home/ysu27/Dropbox/';
    MEGAVEST_path = '/data1/home/ysu27/Dropbox/research/';
else
    prop = 'limits';
    op_path = '/Users/yangsu/Dropbox/';
    MEGAVEST_path = '/Users/yangsu/Dropbox/research/';
%     MEGAVEST_path = 'D:\Dropbox\research\';
end

% ts_cc = tlim;
% datevec(tlim/24/3600+init_time)

% prnlist = 29;
% year = '2014';
% doy = '051';
% flag  = 'ccmin';
% ccmin = 0.6;
% dtau = 60;
matfilestruct = dir([MEGAVEST_path,'xcorr_',year,'_',doy,'*.mat'])
for kk = 1:length(prnlist)
    matfilesname = dir([MEGAVEST_path,'xcorr_',year,'_',doy,'_PRN',num2str(prnlist(kk)),'_*.mat']);
    if ~isempty(matfilesname)
        load([MEGAVEST_path,matfilesname.name]);
    end
end



switch flag
    case 'ccmin'
        ylabelstr = '\rho_c';   
        ytick = 1:length(v_ccmin);
        yticklabel = v_ccmin;
        i_dtau = find(v_dtau==dtau);
        i_ccmin = 1:length(v_ccmin);
        titlestr = (['\Delta\tau = ',num2str(v_dtau(i_dtau)),' s']);
    case 'dtau'
        ylabelstr = '\Delta\tau [s]';   
        ytick = 1:length(v_dtau);
        yticklabel = v_dtau;
        i_ccmin = find(v_ccmin==ccmin);
        i_dtau = 1:length(v_dtau);
        i_dtau = fliplr(i_dtau);
        titlestr = (['\rho_c = ',num2str(v_ccmin(i_ccmin))]);
end

%rotation from geographic drift vg:[vge,vgn,vgz] to magnetic field drift vmf:[vipe,vipn,vi6]
%actually a two-step rotation: 
%geographic -> geomagnetic -> magnetic field line
%declination angle delta, taken from Heinselman and Nicolls, 2008
de = 18.5*pi/180;
%inclination (dip) angle I, taken from Heinselman and Nicolls, 2008
in = 77.5*pi/180;
%rotation matrix
R = [cos(de) -sin(de) 0;
    sin(in)*sin(de) cos(de)*sin(in) cos(in)
    -cos(in)*sin(de) -cos(in)*cos(de) sin(in)];

%i: cc_min idx
%k: dtau idx
%j: time interval idx
% MEGA_VEST2{i,k}(:,j) = [1; 2; 3;    4;  5;  6; 7;    8;   9;     10; 11; 12; 13  14]
% MEGA_VEST2{i,k}(:,j) = [ts;te;V;Psi_v;vge;vgn;AR;Psi_a;flag;npoints;vpe;vpn;vap;vgz]

%drift vamp and vang and vchar [3 4 13]
%anistropy [7 8]
%drfit vperpe and vperpn [11 12]
%drfit vge and vgn [5 6]

% %quiver plot
% for i = 2
%         for k = i_dtau
%             jj = 1:size(MEGA_VEST2{i,k},2);
%             max(MEGA_VEST2{i,k}(3,:))            
%             vmag = MEGA_VEST2{i,k}(col(1),jj);
%             vang = MEGA_VEST2{i,k}(col(2),jj)*pi/180;
%             [vx,vy] = pol2cart(vang,vmag);
% %             quiver(MEGA_VEST2{i,k}(1,jj)/24/3600+init_time,i*ones(size(MEGA_VEST2{i,k}(1,jj))),vx,vy);
%             hold on;
%         end
% end
% xlim(ts_cc/24/3600+init_time);  
% datetick('x','HH:MM:SS','keeplimits');

col = [3 4 13 7 8];
% col = [3 4];
% col = [7 8];
% col = [3 4];

figure('DefaultAxesFontSize',14);
set(gcf,'paperpositionmode','auto');
set(gcf,'units','normalized');

%line plot
for i = i_ccmin
    for j = i_dtau       
        if length(i_dtau) == 1
        switch i
            case 1
                lcolor = 'k';
                marker = 'd';
            case 2
                lcolor = 'r';
                marker = '^';
            case 3
                lcolor = 'c';
                marker = 'o';
            case 4
                lcolor = 'b';
                marker = 's';
            case 5
                lcolor = 'm';
                marker = '.';
            case 6
                lcolor = 'g';
                marker = 'v';
        end
        elseif length(i_ccmin) == 1
        switch j
            case 1
                lcolor = 'k';
                marker = 'd';
            case 2
                lcolor = 'r';
                marker = '^';
            case 3
                lcolor = 'c';
                marker = 'o';
            case 4
                lcolor = 'b';
                marker = 's';
            case 5
                lcolor = 'm';
                marker = '.';
            case 6
                lcolor = 'g';
                marker = 'v';
        end
        end  
        
        for subi = 1:length(col)
        sp(subi) = subplot(length(col),1,subi); 
        switch col(subi)
            case 3
                lblstr = ['(a) SAGA Horizontal Drift Magnitude v_{mag} [m/s], '];
%                 set(sp(subi),'ytick',0:500:1500);
                ylim([0 1500]);
            case 4
                lblstr = ['(b) SAGA Horizontal Drift Orientation v_{ang} [\circ], '];
                set(sp(subi),'ytick',-180:90:180);
%                             set(sp(subi),'yticklabel',{'West','South','East','North','West'});
                ylim([-180 180]);
            case 7
                lblstr = ['(d) SAGA Axial Ratio of Anisotropy Ellipse, '];
            case 8
                lblstr = ['(e) SAGA Anisotropy Ellipse Orientation, \Psi_{\alpha} [\circ], '];
                set(sp(subi),'ytick',0:90:180);
                ylim([0 180]);
            case 5
                lblstr = ['(a) SAGA V_{ge} [m/s], '];
            case 6
                lblstr = ['(b) SAGA V_{gn} [m/s], '];
            case 13
                lblstr = ['(c) SAGA Characteristic Velocity V_{char} [m/s], '];
%                 set(sp(subi),'ytick',0:500:1500);
                ylim([0 1000]);
        end 
        title([lblstr,titlestr]);
        vest = MEGA_VEST{i,j}([1:2 col(subi)],:);
        if col(subi) == 8
            minus_ind = find(vest(3,:)<0);
            vest(end,minus_ind) = 180 + vest(end,minus_ind);
        end
        estimates = vest(end,:);
        estimates = estimates(~isnan(estimates));
        
        dvec = datevec(init_time+tlim/24/3600);
        dvecstt = dvec(1,2:end);
        dvecend = dvec(end,2:end);
        tttt = datestr(init_time+t/24/3600,'HH:MM:SS');
        veststats(subi,1) = mean(estimates);
        veststats(subi,2) = std(estimates);
        
        save([MEGAVEST_path,year,'_',doy,'_'...
            'PRN',num2str(prn),'_estimates','.mat'],...
        'veststats','prn','year','doy','dvecstt','dvecend','tttt');
    
        tc = datenum(mean(vest(1:2,:))/24/3600+init_time);
        time = datenum((vest(1:2,:))/24/3600+init_time);
        for iii = 1:size(vest,2)
            TIME(3*(iii-1)+1:3*iii,:) = [time(1,iii);time(2,iii);NaN]; 
            VEST(3*(iii-1)+1:3*iii,:) = [vest(3,iii);vest(3,iii);NaN];
        end
        if length(i_dtau) == 1
        h_(i) = plot(tc,vest(3,:),[marker,''],'LineWidth',1,'Markersize',6,'color',lcolor,'markerfacecolor',lcolor);
%         h_(i) = plot(TIME,VEST,'LineWidth',1,'color',lcolor); 
        elseif length(i_ccmin) == 1
%         h_(j) = plot(tc,vest(3,:),[marker,''],'LineWidth',1,'Markersize',3,'color',lcolor,'markerfacecolor',lcolor);
        h_(j) = plot(TIME,VEST,'LineWidth',2,'color',lcolor); 
        end
        hold on;grid on;
        
        clear TIME VEST;
        datetick('x');                         
        
        if subi ~= length(col)
            set(sp(subi),'XTickLabel',[]);
        else
        end
        end
    end
end

% t_pfisr = datevec(get(sp(subi),'xlim'));
% % t_pfisr = datevec(tlim/24/3600+init_time);
% flag = 'lat';
% % flag = 'dtau';
% [megadata,lat,dtau] = plotPFISRvs(t_pfisr(1,:),t_pfisr(2,:),flag);
% switch flag
%     case 'dtau'
%         for jj = 1:length(dtau)
%             latdata = megadata(:,1,jj);
%             for jjj = 1:2
%             eb = errorbar(sp(jjj),latdata{jjj}(:,1),latdata{jjj}(:,2),latdata{jjj}(:,3),'x','color','g','Markersize',9,'LineWidth',1);
%             end
%             lgstr = [num2str(lat),'\circ'] ;
%         end
%     case 'lat'
%         for jj = 1:length(lat)
%             latdata = megadata(:,jj,1);
%             for jjj = 1:2
%             eb = errorbar(sp(jjj),latdata{jjj}(:,1),latdata{jjj}(:,2),latdata{jjj}(:,3),'x','color','g','Markersize',9,'LineWidth',1);
%             end
%             lgstr = [num2str(dtau),' s'] ;
%         end
% end

xlabel(['Time [HH:MM UT] on: ',datestr(init_time,'mm/dd/yyyy')]);
% lg = legend(h_,num2str(yticklabel'),'location','North','orientation','horizontal');
if strcmp(flag,'ccmin')
    lg = legend(sp(2),num2str(yticklabel'),'location','North','orientation','horizontal');
else
    lg = legend(sp(2),num2str(yticklabel'),'location','North','orientation','horizontal');
end
% set(findall(gcf,'-property','FontSize'),'FontSize',14);

% for iii = length(col):-1:1
%     op{iii} = get(sp(iii),'outerposition');
%     op{iii}
% end
% %left bottom width height
% set(sp(5),'outerposition',[0    0    1.0000    0.24])
% set(sp(4),'outerposition',[0    0.19    1.0000    0.2])
% set(sp(3),'outerposition',[0    0.34    1.0000    0.2])
% set(sp(2),'outerposition',[0    0.49    1.0000    0.23])
% set(sp(1),'outerposition',[0    0.67   1.0000    0.25])

% get(lg,'Position')
% set(lg,'Position',[0.25 0.93 0.5 0.025]);

return;

figure;
%color-scaled block plot
for subi = 1:length(col)
%     sp(subi) = subplot(2,3,subi);
    sp(subi) = subplot(length(col),1,subi);
    cmap = colormap(jet);
    ylabel(ylabelstr);
    if col(subi) == 4 || col(subi) == 8 || col(subi) == 16
        clim = [-180 180];
    elseif col(subi) == 7
        clim = [0 10];
    else
        clim = [-1500 1500];
    end
    caxis(clim);    
    switch col(subi)
        case 5
            xlabel(['(a) SAGA V_{ge} [m/s], ',titlestr]);
            cblabel = ' drift estimates in geographic coordinates [m/s]';
        case 6
            xlabel(['(b) SAGA V_{gn} [m/s], ',titlestr]);
            cblabel = 'SAGA drift estimates in geographic coordinates [m/s]';
        case 11
            xlabel(['(a) SAGAV_{\perpe} [m/s], ',titlestr]);
            cblabel = 'SAGA drift estimates in magnetic field line coordinates [m/s]';
        case 12
            xlabel(['(b) SAGA V_{\perpn} [m/s], ',titlestr]);
            cblabel = 'SAGA drift estimates in magnetic field line coordinates [m/s]';
        case 3
            xlabel(['(a) SAGA V_{g} [m/s], ',titlestr]);
            cblabel = 'drift [m/s]';
        case 4
            xlabel(['(b) SAGA \Psi_{v_g} [\circ], ',titlestr]);
            cblabel = 'drift [\circ]';
        case 7
            xlabel(['(c) SAGA Axial Ratio, ',titlestr]);
            cblabel = 'anisotropy';
        case 8
            xlabel(['(d) SAGA \Psi_a [\circ], ',titlestr]);
            cblabel = 'anisotropy [\circ]';
        case 15
            xlabel(['(a) SAGA V_{mf} [m/s], ',titlestr]);
            cblabel = 'SAGA drift [m/s]';
        case 16
            xlabel(['(b) SAGA \Psi_{v_{mf}} [\circ], ',titlestr]);
            cblabel = 'SAGA drift [\circ]';
    end
 
    for i = i_ccmin
        for k = i_dtau
            for j = 1:size(MEGA_VEST2{i,k},2)
                xdata = [MEGA_VEST2{i,k}(1,j);MEGA_VEST2{i,k}(2,j);...
                    MEGA_VEST2{i,k}(2,j);MEGA_VEST2{i,k}(1,j);];
                xdata = datenum(xdata/24/3600+init_time);
                if strcmp(flag,'ccmin') == 1
                    %v_ccmin as ytick,v_dtau as ztick
                    ydata = [i;i;i+1;i+1];
                    zdata = ones(4,1)*v_dtau(k);
                elseif strcmp(flag,'dtau') == 1
                    %v_ccmin as ztick,v_dtau as ytick
                    ydata = [k;k;k+1;k+1];
                    zdata = ones(4,1)*v_ccmin(i);
                end
                if MEGA_VEST2{i,k}(9,j)==2 && col(subi)~=13
                    p = patch(xdata,ydata,zdata,'b');
                    set(p,'EdgeColor','w','FaceColor',[0.5 0.5 0.5]);
                elseif MEGA_VEST2{i,k}(9,j)==1 && col(subi)~=13
                    p = patch(xdata,ydata,zdata,'b');
                    set(p,'EdgeColor','w','FaceColor','k');
                else
                    p = patch(xdata,ydata,zdata,'b');
                    set(p,'EdgeColor','w','FaceColor','flat','CData',MEGA_VEST2{i,k}(col(subi),j),'CDataMapping','scaled');
                end
                hold on;   
            end
        end
    end 
    grid on;
    box on;
    set(sp(subi),'YTick',ytick+1,'YTickLabel',yticklabel,'layer','top');
    axis([ts_cc/24/3600+init_time 1 ytick(end)+1]);  
	datetick('x','HH:MM','keeplimits');    
    switch col(subi)
        case {3,4,7,8,15,16}
            cb = colorbar('location','EastOutside',prop,clim);
            set(get(cb,'yLabel'),'String',cblabel);
        case {6,12}
            axpos = get(sp(subi),'position');
            cb = colorbar('location','NorthOutside',prop,clim);
            set(get(cb,'xLabel'),'String',cblabel);
            cbpos = get(cb,'Position');
            set(sp(subi),'Position',[axpos(1) axpos(2)*1.5 axpos(3:end)]);
            set(cb,'Position',[axpos(1) 0.02 axpos(3) 0.02]);
    end
    if subi ~= length(col)
        set(sp(subi),'XTickLabel',[]);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',14);
SAGA_sp = allchild(sp);
SAGA_sp = get(sp,'children');
SAGA_sp = SAGA_sp(1:2,:);

