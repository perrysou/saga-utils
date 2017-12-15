function [peak,tau0,t0] = dataxcorr_alt(sitenum_op,xdata,combos,phstr)

if nargin == 0
    load testdata1.mat
    phstr = 'phase';
    sitenum_op = rx2site(rcvr_op);
    combos = nchoosek(1:size(rcvr_op,1), 2);
end

% specify which phase data are used
if strcmp(phstr,'phase') == 1
    colnum = 3;
    dt = 0.01;
elseif strcmp(phstr,'JROpow') == 1
    colnum = 2;
    dt = 0.02;
elseif strcmp(phstr,'JROph') == 1
    colnum = 3;
    dt = 0.02;
end

% if nargin == 0
%     combos = combos([1 5:7],:);
%     combos(1,:) = fliplr(combos(1,:));
% end

for i = 1:size(combos,1)
    if colnum == 2
        xdata{combos(i,1)}(:,colnum) = detrend(xdata{combos(i,1)}(:,colnum));
        xdata{combos(i,2)}(:,colnum) = detrend(xdata{combos(i,2)}(:,colnum));
    end
    [cc(:,i),~] = xcorr(xdata{combos(i,1)}(:,colnum), xdata{combos(i,2)}(:,colnum));
    [cc2(:,i),lag] = xcorr(xdata{combos(i,2)}(:,colnum), xdata{combos(i,1)}(:,colnum));    
    ac(:,i) = xcorr(xdata{combos(i,1)}(:,colnum));
    ac2(:,i) = xcorr(xdata{combos(i,2)}(:,colnum));
        
    %Normalization   
    cc(:,i) = cc(:,i)/sqrt(max(ac(:,i))*max(ac2(:,i)));
    cc2(:,i) = cc2(:,i)/sqrt(max(ac(:,i))*max(ac2(:,i)));
    ac(:,i) = ac(:,i)/sqrt(max(ac(:,i))*max(ac2(:,i)));
    ac2(:,i) = ac2(:,i)/sqrt(max(ac(:,i))*max(ac2(:,i)));
    
    peak(i) = max(cc(:,i));
    altpeak(i) = cc(lag==0,i);
    altpeak2(i) = cc2(lag==0,i);
    if ~isempty(lag(cc(:,i)==max(cc(:,i)))*dt)
        tau0(i) = lag(cc(:,i)==max(cc(:,i)))*dt;
    else
        continue;
    end
    t0_ind = find(abs(ac(:,i)-max(cc(:,i)))==min(abs(ac(:,i)-max(cc(:,i)))));
%     err_t0(i) = unique(abs(ac(t0_ind)-max(cc(:,i))));
    t0_indm(:,i) = t0_ind(1);
    t0(i) = unique(abs(lag(t0_ind)*dt));  
end

[~,order_t] = sort(tau0,'ascend');
% [~,order_cc] = sort(peak,'descend');
order = order_t;

cmap = hsv(size(combos,1));

for i = 1:size(combos,1)
    j = order(i);
%     text(min(lag)*dt*0.95,1+(size(combos,1)-i)*2,...
%         {['$\rho_{max}=$',num2str(peak(j),'%.4f')];...
% %         ['@$\tau_{cm}=$',num2str(tau0(j),'%.2f'),'s,'];...
% %         ['$\tau_{am}=$',num2str(t0(j),'%.2f'),'s']},...
%         'VerticalAlignment','Top');
% %         '$\rho(0)=$',num2str(altpeak(j),'%.4f')],...
% %         'VerticalAlignment','Top');
    hold on;
    h(i) = plot(lag*dt, cc(:,j)+(size(combos,1)-i)*2,'Color',cmap(i,:)); 
%     drawnow;
    plot(lag*dt, ac(:,j)+(size(combos,1)-i)*2,'Color','c'); 
    if tau0(j) <= 0
        plot([tau0(j);t0(j)],[peak(j);ac(t0_indm(j),j)]+(size(combos,1)-i)*2,'k');
    else
        plot([tau0(j);-t0(j)],[peak(j);ac(t0_indm(j),j)]+(size(combos,1)-i)*2,'k');    
    end

    lg{i,:} = [sitenum_op{combos(j,1),:},'\&',sitenum_op{combos(j,2),:}];
%     sitenum_op{combos(j,1),:},'&',sitenum_op{combos(j,2),:}
%     plot(0,altpeak(j)++(i-1)*2,'k.');
end 
legend(h, lg);
set(gca,'Layer','top','XGrid','on','YGrid','on','Ytick',-1:2:(size(combos,1)-1)*2+1);
set(gca,'YtickLabel',[]);
ylim([-1 (size(combos,1)-1)*2+1]);
ylabel('Normalized Cross-correlation Coefficients');
ylabel('Normalized $\rho_{ij}$');
xlabel('Lag time [s]');
title(['$\rho_{max}=$',num2str(peak(j),'%.4f'),...
        '@$\tau_{cm}=$',num2str(tau0(j),'%.2f'),'s,',...
        '$\tau_{am}=$',num2str(t0(j),'%.2f'),'s']);
% title('Cross-correlation over 2615-2660s after 0300UT for IIT-3');
% [~,op_path] = ver_chk;
% plotpath = [op_path, 'PRN23_Lag_plot_2615-2660s_after_0300UT_iit3', '.eps'];
% saveas(gcf,plotpath,'epsc2');
end

