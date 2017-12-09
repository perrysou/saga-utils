function [] = diststats(lrfile)
%estimate the statistical distribution of the scintillation data
close all;
load 'sd_2014.mat';
MSP = sortrows(MEGA_MSP,1);

% MSP = MS4;
MSP = sortrows(MSP,2);
% plot(MSP(:,2));
max(MSP(:,2))
[phat_r,pci_r] = mle(MSP(:,2),'distribution','rayleigh')
[phat_n,pci_n] = mle(MSP(:,2),'distribution','nakagami')
[phat_g,pci_g] = mle(MSP(:,2),'distribution','normal')
[phat_ga,pci_ga] = mle(MSP(:,2),'distribution','gamma','alpha',0.10)

pdf_r = pdf('rayleigh',MSP(:,2),phat_r);
pdf_n = pdf('nakagami',MSP(:,2),phat_n(1),phat_n(2));
pdf_g = pdf('normal',MSP(:,2),phat_g(1),phat_g(2));
pdf_ga = pdf('gamma',MSP(:,2),phat_ga(1,1),phat_ga(1,2));

[f,x] = ecdf(MSP(:,2));
% stairs(x,f);
pdf_ = ecdfhist(f,x,linspace(0,max(MSP(:,2)),100));
stairs(linspace(0,max(MSP(:,2)),100),pdf_,'k');
% set(get(gca,'Children'),'FaceColor','w')
hold on;
plot(MSP(:,2),pdf_r,'-b');
plot(MSP(:,2),pdf_n,'-c');
plot(MSP(:,2),pdf_g,'-r');
plot(MSP(:,2),pdf_ga,'-g');
legend('actual data','rayleigh','nakagami','normal','gamma');
end

