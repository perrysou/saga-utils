function [] = plot_STEC()
prn = 31;
doy = '342';
year = '2013';
load ionodata.mat;
if ~isempty(IONODATA)
%     t_iono = (datenum(gps2utc(IONODATA(:,1:2)))-init_time)*24*3600;
%     iii = t_iono<=newtelist(tt) & t_iono>=newtslist(tt) & IONODATA(:,5)==prn;
%     STEC(:,rr) = mean(IONODATA(iii,3));   
%     unique(IONODATA(:,5))
    ind = IONODATA(:,5)==prn;
    t_iono = datenum(gps2utc(IONODATA(ind,1:2)));
    STEC = IONODATA(ind,3);
    plot(t_iono,STEC);
end
datetick('x','HH:MM','keepticks')
% legend(sitenum_op);
title(['Slant TEC for PRN:',num2str(prn),' on DOY:',doy,', year:',year]);
% xlabel(['Time after ',datestr(init_time,'HH:MM'),...
% ' UT ',datestr(init_time,'mm/dd/yyyy')]);  
% plotpath = [op_path,'PRN',num2str(prn),'_STEC','.eps'];
% saveas(gcf,plotpath,'epsc2');
% close;
end