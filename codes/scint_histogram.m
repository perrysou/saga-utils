function [] = scint_histogram(lrfile) 
op_path = '/Users/yangsu/Dropbox/';
% lrfile = [op_path,'lrdata_',year,'_',doy,'.mat'];
lrfile = 'lrdata_2013_342.mat'
load(lrfile);
tlim_vec = datevec(lrdata{5})
MSP = lrdata{1};HITDATA = lrdata{2};
CST = lrdata{3};rcvr_op = lrdata{4};
tlim = lrdata{5};splim = lrdata{6};
signal = lrdata{7};spmask = lrdata{8};
MS4 = lrdata{9};s4mask = lrdata{10};
s4lim = lrdata{11};
span = 60;
scintdata = [MSP(:,1) MS4(:,2) MSP(:,2:end)];
size(scintdata)
time = scintdata(:,1);
s4 = scintdata(:,2);
max(s4)
sp = scintdata(:,3);
max(sp)
prn = unique(scintdata(:,4))';
tdata0 = scintdata(find(sp<0.6),1);
tdata1 = scintdata(find(s4<0.2&sp>=0.6),1);
tdata2 = scintdata(find(s4>=0.2&sp>=0.6),1);
% histogram(tdata0,'normalization','count','binwidth',1/24);
hold on;
histogram(tdata1,'normalization','count','binwidth',1/24,'facecolor','w');
hold on;
histogram(tdata2,'normalization','count','binwidth',1/24,'facecolor','k');
datetick('x','HH');
set(gca,'layer','top',...
    'yticklabel',num2str([get(gca,'ytick')*100/size(scintdata,1)]','%.2f'))
xlabel('Hour [UT]');
ylabel('Occurence [%]');
title({'Amplitude S_4 and Phase \sigma_{\phi} Scintillation Events on 12/08/2013',...
    ['S_{4_{max}} = ',num2str(max(s4)),', \sigma_{{\phi}_{max}} = ',num2str(max(sp)),' rad']});
grid on;
legend('S_4<0.2, \sigma_{\phi}>=0.6','S_4>=0.2, \sigma_{\phi}>=0.6');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
% histogram(scintdata(:,2))
% [hdata, grid] = hist3(scintdata(:,2:3));
% [hdata, grid] = hist3(scintdata(:,2:3),'edges',{[0 0.2],[0 0.6 2*pi]});
% hist3(scintdata(:,2:3),'edges',{[0 0.2],[0 0.6]});
% xb = linspace(min(scintdata(:,2)),max(scintdata(:,3)),size(hdata,1));
% yb = linspace(min(scintdata(:,3)),max(scintdata(:,3)),size(hdata,1));
% [X,Y] = meshgrid(grid{1},grid{2});
% contourf(X,Y,hdata,'showtext','on')
% colormap(jet)
% set(get(gca,'child'),'Facecolor','interp','cdatamode','auto')
% set(gca,'zscale','linear')
% scintdata(find(s4>0.2&sp>0.6),2)
% [mega_t,TSP_hr0,TSP_hrv0] = find_general_times(MS4,rcvr_op,0.2);
% TSP_hrv0
% [mega_t,TSP_hr0,TSP_hrv0] = find_general_times(MSP,rcvr_op,0.6);
% TSP_hrv0
saveas(gca,[op_path,'histogram.png'],'png');
% close;
end