function [AZb,ELb,beamid]= plotPFISR_NeTe(t0,tf,flag)
% t0 = [2015 3 17 13 31 24];
% tf = [2015 3 17 13 33 32];
flag = 'Ne';
load('PFISR_data2.txt');

PFISR_data = PFISR_data2;data = [datenum(PFISR_data(:,1:6)) PFISR_data(:,7:end)];
data = data(data(:,1)<=datenum(tf)+10/24/3600&data(:,1)>=datenum(t0)-10/24/3600,:);
if ~isempty(data)
    beamid = unique(data(:,4),'stable');
    for i = 1:length(beamid)
        AZb(i,:) = unique(data(data(:,4)==beamid(i),2));
        ELb(i,:) = unique(data(data(:,4)==beamid(i),3));
    end
else
    beamid = NaN;
    AZb = NaN;
    ELb = NaN;
end

end