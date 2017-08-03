function [ ] = plot_vel(VEST,init_time,ts_proced,lagmin,ccmin_easy,ccmin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% load easyvest.mat;

% %time in datenum
% time = init_time+VEST(1:2,:)/24/3600;
%time in seconds after init_time
time = VEST(1:2,:);
time_u = VEST(1:2,VEST(9,:)==2);
time_i = VEST(1:2,VEST(9,:)==1);

tmid = cumsum(time,1)/2;
tmid = tmid(2,:)';
% %tlim = datenum
% tlim = init_time+ts_proced/24/3600;
% tlim in seconds after init_time
tlim = ts_proced;

for jj = 1:size(time,2)
    TIME(3*(jj-1)+1:3*jj,:) = [time(1,jj)';time(2,jj)';NaN];
end

if size(time_u,2) ~= 0
    for jj = 1:size(time_u,2)
        TIME_U(3*(jj-1)+1:3*jj,:) = [time_u(1,jj)';time_u(2,jj)';NaN];
    end
else
    TIME_U = zeros(size(time_u,2),1);
end
DATA_U = zeros(size(TIME_U));

if size(time_i,2) ~= 0
    for jj = 1:size(time_i,2)
        TIME_I(3*(jj-1)+1:3*jj,:) = [time_i(1,jj)';time_i(2,jj)';NaN];    
    end
else
    TIME_I = zeros(size(time_i,2),1);
end
DATA_I = zeros(size(TIME_I));

ve = VEST(3,:)';
vn = VEST(4,:)';
vec = VEST(5,:)';
vnc = VEST(6,:)';
ar = VEST(7,:)';
psi_a = VEST(8,:)';
for ii = 1:length(ve)
    VE(3*(ii-1)+1:3*ii,:) = [ve(ii);ve(ii);NaN];
    VN(3*(ii-1)+1:3*ii,:) = [vn(ii);vn(ii);NaN];
    VEC(3*(ii-1)+1:3*ii,:) = [vec(ii);vec(ii);NaN];
    VNC(3*(ii-1)+1:3*ii,:) = [vnc(ii);vnc(ii);NaN];
    AR(3*(ii-1)+1:3*ii,:) = [ar(ii);ar(ii);NaN];
    PSI_A(3*(ii-1)+1:3*ii,:) = [psi_a(ii);psi_a(ii);NaN];
end

%time intervals with no available results
tmid_u = tmid(isnan(vec)&VEST(9,:)'==2);
%time intervals with invalid results (a*b-h^2<=0)
tmid_i = tmid(isnan(vec)&VEST(9,:)'==1);
data_u = zeros(size(tmid_u));
data_i = zeros(size(tmid_i));

subplot(4,1,1);
x1 = plot(TIME,VE,'-k',TIME,VEC,'-k',TIME_U,DATA_U,'-k',TIME_I,DATA_I,'-k');
hold on;
y1 = plot(tmid,ve,'.-b',tmid,vec,'d-r',tmid_u,data_u,'sg',tmid_i,data_i,'om');
legend(y1,{'\vee_e','\vee_{c_v}','\vee_{c_u}','\vee_{c_i}'},...
    'Location','NorthEastOutside','Orientation','Vertical');
grid on;
xlim(tlim);
set(gca,'XMinorTick','on');
if min(VE)~=max(VE) && ~isnan(min(VE))
    ylim([min(VE)*1.1 max(VE)*1.1]);
end
% datetick('x','HH:MM','keeplimits');
ylabel('V_{East} (m/s)');
title({'Estimated 2D Drift Velocity';...
    ['abs(tpeak_{min})=',num2str(lagmin),...
    ', ccpeak_{min}=',num2str(ccmin_easy),...
    ', cc_{min}=',num2str(ccmin)]});

subplot(4,1,2);
x2 = plot(TIME,VN,'-k',TIME,VNC,'-k',TIME_U,DATA_U,'-k',TIME_I,DATA_I,'-k');
hold on;
y2 = plot(tmid,vn,'.-b',tmid,vnc,'d-r',tmid_u,data_u,'sg',tmid_i,data_i,'om');
legend(y2,{'\vee_e','\vee_{c_v}','\vee_{c_u}','\vee_{c_i}'},...
    'Location','NorthEastOutside','Orientation','Vertical');
grid on;
xlim(tlim);
set(gca,'XMinorTick','on');
if min(VN)~=max(VN) && ~isnan(min(VN))
    ylim([min(VN)*1.1 max(VN)*1.1]);
end
ylabel('V_{North} (m/s)');

subplot(4,1,3);
x3 = plot(TIME,AR,'-k',TIME_U,DATA_U,TIME_I,DATA_I,'-k');
hold on;
y3 = plot(tmid,ar,'d-r',tmid_u,data_u,'sg',tmid_i,data_i,'om');
ylabel('Axial Ratio');
grid on;
xlim(tlim);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
legend(y3,{'AR_v','AR_u','AR_i'},...
    'Location','NorthEastOutside','Orientation','Vertical');

subplot(4,1,4);
x4 = plot(TIME,PSI_A,'-k',TIME_U,DATA_U,TIME_I,DATA_I,'-k');
hold on;
y4 = plot(tmid,psi_a,'d-r',tmid_u,data_u,'sg',tmid_i,data_i,'om');
ylabel('\Psi_a [deg]');
grid on;
% datetick('x','HH:MM','keeplimits');
xlim(tlim);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
xlabel(['Time after ',datestr(init_time,'HH:MM'),...
    ' UT ',datestr(init_time,'mm/dd/yyyy')]);
legend(y4,{'\Psi_{a_v}','\Psi_{a_u}','\Psi_{a_i}'},...
    'Location','NorthEastOutside','Orientation','Vertical');
end

