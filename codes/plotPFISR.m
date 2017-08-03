function [PFISR_sp,lat,dtau] = plotPFISR(ts,te,flag)
% 
ts = [2014 2 20 11 20 0];
te = [2014 2 20 11 57 0];
flag = 'lat';

if verLessThan('matlab','8.4.0')
    prop = 'clim';
    op_path = '/data1/home/ysu27/Dropbox/';
else
    prop = 'limits';
    op_path = '/Users/yangsu/Dropbox/';
%     op_path = 'D:\Dropbox\';
end

load('PFISR_data1.txt');

PFISR_data = PFISR_data1;
data = [datenum(PFISR_data(:,1:6)) PFISR_data(:,7:end)];
data(:,11) = -data(:,11)+90;
idx = find(data(:,11)>180);
data(idx,11) = data(idx,11)-360;
max(data,11);
data = data(data(:,1)<=datenum(te)+10/24/3600&data(:,1)>=datenum(ts)-10/24/3600,:);
tdt = unique(data(:,1),'stable');
dtau = unique(data(:,2),'rows');
lat = unique(data(:,3),'rows');

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

%inverse transform to geographic velocities
vg = R\data(:,4:6)';
vgmag = sqrt(vg(1,:).^2+vg(2,:).^2);
vgang = 180/pi*atan2(vg(2,:),vg(1,:));
data = [data vg' vgmag' vgang'];
%geographic velocities vge vgn vgz vgmag vgang
%[14 15 16 17 18]
switch flag
    case 'lat'
        ylabelstr = 'Latitude [\circ]';   
        ytick = 1:length(lat);
        yticklabel = lat;
        dtau = dtau(1);        
        titlestr = (['\Delta\tau = ',num2str(dtau),' s']);
    case 'dtau'
        ylabelstr = '\Delta\tau [s]';   
        ytick = 1:length(dtau);
        yticklabel = dtau;
        lat = 66.25;
        titlestr = (['Latitude = ',num2str(lat),'\circ']);
end

%plot vamp and vpsi [10 11]
% %plot vperpe and vperpn [4 5]
% %plot vge and vgn [14 15]
col = [17 18];
errcol = [12 13];
% col = [12 13];
%quiver plot
% for jj = 1:length(dtau)
%     for ii = 1:length(tdt)
%         for kk = 1:length(lat)
%             block = data(data(:,1)==tdt(ii)&data(:,2)==dtau(jj)&data(:,3)==lat(kk),:);
%             if isempty(block)
%                 vmag(kk,ii) = NaN;
%                 vang(kk,ii) = NaN;
%             else
%                 vmag(kk,ii) = block(:,col(1));
%                 vang(kk,ii) = (-block(:,col(2))+90)*pi/180;
%             end
%         end
%     end
%     [vx,vy] = pol2cart(vang,vmag);
%     [tgrid,latgrid] = meshgrid(tdt,lat);
% %     quiver(tgrid,latgrid,vx,vy);
%     datetick('x','HH:MM:SS');
%     hold on;
% end

for kk = 1:length(lat)
    for dt = 1:length(dtau)
        switch dt
        case 1
            linecolor = 'r';
        case 2
            linecolor = 'k';
        case 3
            linecolor = 'cyan';
        case 4
            linecolor = 'r';
        case 5
            linecolor = 'g';
        case 6
            linecolor = 'magenta';   
        otherwise   
            linecolor = 'y';
        end
        for subi = 1:2
            sp(subi) = subplot(2,1,subi);
            latdata = data(data(:,3)==lat(kk)&data(:,2)==dtau(dt),[1 col(subi) errcol(subi)]);
%             plot(latdata(:,1),latdata(:,2),'.-','LineWidth',2,'Markersize',9,'color',linecolor);
            errorbar(latdata(:,1),latdata(:,2),latdata(:,3),'.','color',linecolor);
            switch col(subi)
                case {10,17}
                    title(['(a) PFISR V [m/s], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates [m/s]';
                    set(sp(subi),'ytick',0:1000:5000);
                    ylim([0 5000]);
                case {11,18}
                    title(['(b) PFISR \Psi_v [\circ], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates [\circ]';
                    set(sp(subi),'ytick',-180:90:180);
                    ylim([-180 180]);
                case 4
                    title(['(a) PFISR V_{\perpe} [m/s], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates in magnetic field line coordinates [m/s]';
                    set(sp(subi),'ytick',-5000:1000:5000);
                    ylim([-5000 5000]);
                case 14
                    title(['(a) PFISR V_{ge} [m/s], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates in geographic coordinates [m/s]';
                    set(sp(subi),'ytick',-5000:1000:5000);
                    ylim([-5000 5000]);
                case 5
                    title(['(b) PFISR V_{\perpn} [m/s], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates in magnetic field line coordinates [m/s]';
                    set(sp(subi),'ytick',-5000:1000:5000);
                    ylim([-5000 5000]);
                case 15
                    title(['(b) PFISR V_{gn} [m/s], ',['Latitude = ',num2str(lat(kk)),'\circ']]);
                    cblabel = 'PFISR drift estimates in geographic coordinates [m/s]';
                    set(sp(subi),'ytick',-5000:1000:5000);
                    ylim([-5000 5000]);
            end
            hold on;
            grid on;
            box on;
            xlim([datenum(ts) datenum(te)]);
            datetick('x','HH:MM','keeplimits');  
            set(sp(subi),'layer','top'); 
            xlabel('Universal Time [HH:MM]');
        end
    end
lg = legend(cellstr(num2str(dtau)),'Orientation','Horizontal');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
get(lg,'Position');
set(lg,'Position',[0.4 0.01 0.2 0.05]);
saveas(gcf,[op_path,'PFISR_vest_',num2str(lat(kk)),'s_'...
    'errorbar',...
    '.eps'],'epsc2');
  close;
end
return;

figure;
col = [col errcol];
% keyboard;
for subi = 1:length(col)
    sp(subi) = subplot(length(col)/2,2,subi);
    colormap(jet);
    if col(subi) == 11 || col(subi) == 18 || col(subi) == 13
        clim = [-180 180];
    else
        clim = [-1500 1500];
    end
    caxis(clim);
    
    switch col(subi)
    case 10
        title(['(a) V_{mf} [m/s], ',titlestr]);
        cblabel = 'PFISR drift estimates [m/s]';
    case 11
        title(['(b) \Psi_{v_{mf}} [\circ], ',titlestr]);
        cblabel = 'PFISR drift estimates [\circ]';
    case 17
        title(['(a) V_{g} [m/s], ',titlestr]);
        cblabel = 'PFISR drift [m/s]';
    case 18
        title(['(b) \Psi_{v_{g}} [\circ], ',titlestr]);
        cblabel = 'PFISR drift [\circ]';
    case 4
        title(['(a) V_{\perpe} [m/s], ',titlestr]);
        cblabel = 'PFISR drift estimates in magnetic field line coordinates [m/s]';
    case 5
        title(['(b) V_{\perpn} [m/s], ',titlestr]);
        cblabel = 'PFISR drift estimates in magnetic field line coordinates [m/s]';
    case 12
        title(['(c) eV_{mf} [m/s], ',titlestr]);
        cblabel = 'PFISR error [m/s]';
    case 13
        title(['(d) e\Psi_{v_{mf}} [\circ], ',titlestr]);
        cblabel = 'PFISR error [\circ]';
    end
    
    for ii = 1:length(tdt)
        for jj  = 1:length(dtau)
            for kk = 1:length(lat)
                if tdt(ii)+dtau(jj)/24/3600 <= datenum(te)
                    xdata = [tdt(ii);tdt(ii)+dtau(jj)/24/3600;tdt(ii)+dtau(jj)/24/3600;tdt(ii)];
                else
                    xdata = [tdt(ii);datenum(te);datenum(te);tdt(ii)];                    
                end
                switch flag
                    case 'dtau'
                    ydata = [jj;jj;jj+1;jj+1];
                    zdata = ones(4,1)*lat(kk);
                    case 'lat'
                    ydata = [kk;kk;kk+1;kk+1];
                    zdata = ones(4,1)*dtau(jj);  
                end
                data_block = data(data(:,1)==tdt(ii)&data(:,2)==dtau(jj)&data(:,3)==lat(kk),col(subi));
                if ~isempty(data_block)
                    p = patch(xdata,ydata,zdata,'b');
                    set(p,'EdgeColor','w','FaceColor','flat','CData',data_block,'CDataMapping','scaled');
                end
                hold on;
            end
        end
    end
    grid on;
    box on;
    ylabel(ylabelstr);
    set(sp(subi),'ytick',ytick+1,'yticklabel',yticklabel,'layer','top');
    if subi == 2 || subi == 4
        set(sp(subi),'YTickLabel',[]);
        ylabel([]);
    end       
    xlim([datenum(ts) datenum(te)]);
    ylim([ytick(1) ytick(end)+1]);
    datetick('x','HH:MM','keeplimits');
    switch col(subi)
        case {10,11,14,15,17,18,12,13}
            cb = colorbar('location','EastOutside',prop,clim);
            set(get(cb,'yLabel'),'String',cblabel);
        case {5,6}
            axpos = get(sp(subi),'position');
            cb = colorbar('location','NorthOutside',prop,clim);
            set(get(cb,'xLabel'),'String',cblabel);
            cbpos = get(cb,'Position');
            set(sp(subi),'Position',[axpos(1) axpos(2)*1.5 axpos(3:end)]);
            set(cb,'Position',[axpos(1) 0.02 axpos(3) 0.02]);
            set(sp(subi),'XTickLabel',[]);
    end

end

set(findall(gcf,'-property','FontSize'),'FontSize',14);

PFISR_sp = allchild(sp);
PFISR_sp = get(sp,'children');
PFISR_sp = PFISR_sp(1:2,:);
end

