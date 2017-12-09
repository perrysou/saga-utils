function plotmisc(xcorr_results,year,doy,prn)
%initialization
if verLessThan('matlab','8.4.0')
    prop = 'clim';
    op_path = '/data1/home/ysu27/Dropbox/';
    load(xcorr_results);
else
    prop = 'limits';
    op_path = '/Users/yangsu/Dropbox/';
end

% %set time axis limit for PFISR
% t_pfisr = datevec(tlim/24/3600+init_time)

%plot saga vest with fixed dtau
figure;
plot_vel_v3('dtau',[],0.6,year,doy,prn)
plotname = ['PRN',num2str(prn),'_SAGA_vest_',year,'_',doy...
    '_zoom',num2str(zcounter)','_'...
    num2str(tlim(1)),'-',num2str(tlim(2)),'s_after_',...
    datestr(init_time,'HHMM'),'UT_','dtau'];
plotpath = [op_path,plotname,'.eps'];
saveas(gcf,plotpath,'epsc2');
close;

return;

%plot saga vest with fixed ccmin
figure;
plot_vel_v3('ccmin',60,[],year,doy,prn)
plotname = ['PRN',num2str(prn),'_SAGA_vest_',year,'_',doy,...
    '_zoom',num2str(zcounter)','_'...
    num2str(tlim(1)),'-',num2str(tlim(2)),'s_after_',...
    datestr(init_time,'HHMM'),'UT_','ccmin'];
plotpath = [op_path,plotname,'.eps'];
saveas(gcf,plotpath,'epsc2');
close;

return;

figure;
%plot pfisr with fixed dtau
[PFISR_sp,lat,dtau_pfisr] = plotPFISR(t_pfisr(1,:),t_pfisr(2,:),'lat')
plotname = ['PRN',num2str(prn),'_PFISR_vest_',year,'_',doy,...
    '_zoom',num2str(zcounter)','_'...
    num2str(tlim(1)),'-',num2str(tlim(2)),'s_after_',...
    datestr(init_time,'HHMM'),'UT_','lat'];
plotpath = [op_path,plotname,'.eps'];
saveas(gcf,plotpath,'epsc2');

figure;
%plot pfisr with fixed lat
[PFISR_sp,lat,dtau_pfisr] = plotPFISR(t_pfisr(1,:),t_pfisr(2,:),'dtau')
plotname = ['PRN',num2str(prn),'_PFISR_vest_',year,'_',doy,...
    '_zoom',num2str(zcounter)','_'...
    num2str(tlim(1)),'-',num2str(tlim(2)),'s_after_',...
    datestr(init_time,'HHMM'),'UT_','dtau'];
plotpath = [op_path,plotname,'.eps'];
saveas(gcf,plotpath,'epsc2');

%pfisr & saga vmf comparison
figure;
colormap(jet);
lat_str = [', Lat = ',num2str(lat),' \circ'];
dtau_str = [', \Delta\tau = ',num2str(v_dtau(1)),' s'];
for subi = 1:4
    vg(subi) = subplot(4,1,subi);
    grid on;  
    box on;
    set(vg(subi),'layer','top');
    switch subi
        case {1,3} 
            copyobj(PFISR_sp{(subi+1)/2},vg(subi));
            ylabel('\Delta\tau [s]');
            set(vg(subi),'YTick',(1:length(dtau_pfisr))+1,'YTickLabel',dtau_pfisr);           
            axis([datenum(t_pfisr(1,:)) datenum(t_pfisr(2,:)) 1 length(dtau_pfisr)+1]);
        case {2,4}          
            copyobj(SAGA_sp{subi/2},vg(subi));
            ylabel('N_p');
            set(vg(subi),'YTick',(1:length(v_ccmin))+1,'YTickLabel',v_ccmin);
            axis([tlim/24/3600+init_time 1 length(v_ccmin)+1]);
    end
    switch subi
        case 1
            clim = [-1500 1500];
            cblabel = ['(a) PFISR V [m/s]',lat_str];
        case 2
            clim = [-1500 1500];
            cblabel = ['(b) SAGA V [m/s]',dtau_str];
        case 3
            clim = [-180 180];
            cblabel = ['(c) PFISR \Psi_v [\circ]',lat_str];
        case 4
            clim = [-180 180];
            cblabel = ['(d) SAGA \Psi_v [\circ]',dtau_str];
    end
    caxis(clim);
    xlabel(cblabel);
    cb = colorbar('location','EastOutside',prop,clim);
%     set(get(cb,'yLabel'),'String',cblabel);
    datetick('x','HH:MM','keeplimits');
    if subi ~= 4
        set(vg(subi),'XTickLabel',[]);
    end
end

set(findall(gcf,'-property','FontSize'),'FontSize',14);
% set(vmf,'FontSize',12);

plotname = ['PRN',num2str(prn),'_V_g_',...
num2str(tlim(1)),'-',num2str(tlim(2)),'s_after_',...
datestr(init_time,'HHMM'),'UT'];
plotpath = [op_path,plotname,'.eps'];
saveas(gcf,plotpath,'epsc2');
close all;