clear all;
close all;
load sd_2014.mat;
% load sd.mat;
op_path = '/Users/yangsu/Dropbox/research-misc/ION2015_slides/';
%%
% subplot(2,1,1);
truetime = MSP_days(:,1)+MEGA_MSP(:,1)-1;
ttv = datevec(truetime);
%local time
truetime = truetime-9/24;
ttv = datevec(truetime);
id = find(truetime>=datenum([2014 1 1 0 0 0]));
truetime = truetime(id,:);

MEGA_MSP = MEGA_MSP(id,:);
ttv = datevec(truetime);
doy = datenum([ttv(:,1:3) zeros(size(ttv,1),3)])-datenum([2014 1 1 0 0 0])+1;
toffset = truetime-doy+1;
ttv = datevec(toffset);
MSP = [doy toffset MEGA_MSP(:,2:end)];
colormap(jet);

% hold on;
% for i = 50:52
%     MSP_day = MSP(doy==i,:);
%     scatter(MSP_day(:,1)+datenum([2014 1 1 0 0 0])-1,MSP_day(:,2),[],MSP_day(:,3),'.');
%     clim = [0 2*pi];
%     caxis(clim);
% end
% cb = colorbar;
% % get(cb,'limits')
% set(get(cb,'YLabel'),'String','\sigma_\Phi value [Rad]');
% if verLessThan('matlab','8.4.0')
%     set(cb,'Clim',clim);
% else
%     set(cb,'limits',clim);
% end
% datetick('y','HH:MM');
% datetick('x','mm-dd');
% xlabel('day of year 2014');

% subplot(2,1,2);
hold on;
for i = 50:51
    MSP_day = MSP(doy==i,:);
%     scatter(MSP_day(:,1)+datenum([2014 1 1 0 0 0])-1,MSP_day(:,2),[],MSP_day(:,3),'.');
    sc = scatter3(MSP_day(:,1)+datenum([2014 1 1 0 0 0])-1,MSP_day(:,2),MSP_day(:,3),9,MSP_day(:,3),'.');
    clim = [0 2*pi];
    caxis(clim);
end
cb = colorbar;
% get(cb,'limits')
set(get(cb,'YLabel'),'String','\sigma_\Phi value [Rad]');
if verLessThan('matlab','8.4.0')
    set(cb,'Clim',clim);
else
    set(cb,'limits',clim);
end

%plot daily solar terminators
for jjj = 1:365
dstr = datestr(jjj+datenum([2014 1 1 0 0 0])-1,'yyyy-mm-dd');
RH_(jjj,:) = datenum(['2014-01-01 ',RH(jjj,:)],'yyyy-mm-dd HHMM');
SH_(jjj,:) = datenum(['2014-01-01 ',SH(jjj,:)],'yyyy-mm-dd HHMM');
end
rh = plot3([1:365]+datenum([2014 1 1 0 0 0])-1,RH_,7*ones(size(RH_)),'r','Linewidth',2);
sh = plot3([1:365]+datenum([2014 1 1 0 0 0])-1,SH_,7*ones(size(SH_)),'k','Linewidth',2);
hold on;
legend([rh;sh],'Rise time','Set time');
% title('Solar terminator for year 2014');

datetick('y','HH:MM');
datetick('x','mmm');
xlabel('Day of Year 2014');
ylabel('Local Time HH:MM');
box on;grid on;
set(gca,'layer','top');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
% view([-30 30]);
view([-0 90]);
% saveas(gcf,[op_path,'dailySPwsolarterm.jpg'],'jpg');
% close;
