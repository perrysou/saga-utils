function [] = hr_routine(yearin,doyin)
%% Initialization
close all;
    
% maxNumCompThreads(8);

%file separator "/" in linux
sep = filesep;
%specify signal 
signal_type = 0;
signal = 'L1CA';

%where cases data are stored
cases_folder = '/data1/public/Data/cases/pfrr/';
%case_folder = '/data1/public/Data/cases/calg/';

%where output data and plots are stored
home_dir = '/data1/home/ysu27/';

%year
% year = '2014';

%find yesterday in doy
comm = 'date -u -d "a day ago" +%j';
[~,yesterday] = system(comm);
yesterday = str2num(yesterday);

comm1 = 'date -u -d "a day ago" +%Y';
[~,year] = system(comm1);
yearnum = str2num(year);

if nargin == 2
    year = num2str(yearin,'%04i');
    doylist = doyin;
else
    year = num2str(yearnum,'%04i');
    doylist = yesterday;
end

% for doy = yesterday
for doynum = doylist
%day of year in string
doy = num2str(doynum,'%03i');
%% Process low rate data 

% % Low-rate data processing(can be skipped if done before) 
% %sigmaphi mask in [rad]
% spmask = 0.5;
% %generate low rate sigmaphi plots for all operational receivers for the doy
% %with sigmaphi mask
% [MSP,HITDATA,CST,rcvr_op,tlim,splim] = scint_el_stackplot(cases_folder,home_dir,signal_type,doy,year,spmask);
% lrdata = {MSP,HITDATA,CST,rcvr_op,tlim,splim};
if strcmp(cases_folder(end-4:end-1),'pfrr')
    %folder_path for 2013 Poker Flat data 
    op_path = [home_dir,'PFRR_Data/'];
else
    %folder_path for 2013 Calgary data
    op_path = [home_dir,'Calgary_Data/'];  
end
% save([op_path,'lrdata_',year,'_',doy,'.mat'],'lrdata');
%% Plot low-rate sigmaphi stackplots 
%load low rate results
load([op_path,'lrdata_',year,'_',doy,'.mat']);
% load('/media/yang/Misc/Dropbox/research/lrdata_2013_342.mat');
hr_times = [];
hrtimes = [];
MSP = lrdata{1};HITDATA = lrdata{2};
CST = lrdata{3};rcvr_op = lrdata{4};
tlim = lrdata{5};splim = lrdata{6};
span = 60;
%receiver structure including all poker flat receivers
rcvr_struct = rx_dirs(cases_folder,year,doy);

% if ~isempty(MSP)
%     lr_rx_stackplot(lrdata,year,doy,rcvr_struct,op_path);
%     lr_prn_stackplot(lrdata,year,doy,rcvr_struct,op_path);
% end
% 
% %% Find interesting times for high-rate processing
% % PRNs detected to have scintillation
% %unique(CST(:,2))'
% if ~isempty(CST)&&~isempty(HITDATA)
% %     for prn = unique(CST(:,2))'
%     for prn = unique(CST(:,2))'
%         MAX = [];    
%         for rr = unique(HITDATA(:,4))'
%             hitdata = HITDATA(HITDATA(:,4)==rr&HITDATA(:,3)==prn,:);
%             maxdata = hitdata(hitdata(:,2)==max(hitdata(:,2)),:);
%             MAX = [MAX;maxdata];
%         end
%         [tt,dt] = find_tspan(MAX,span); 
%         HHMM = datestr(tt,'HHMM');
% %         hr_times = [hr_times;prn,str2num(year),str2num(doy),str2num(HHMM(1,:)),str2num(HHMM(2,:))];
%         hr_times = [hr_times;num2str(prn,'%02i'),' ',year,' ',doy,' ',HHMM(1,:),' ',HHMM(2,:)];
%         hrtimes = [hrtimes;prn tt' dt];
%     end   
% end


%% Discard scintillation values less than a larger sigmaphi threshold
spth = 0.8;
if ~isempty(MSP)
    if ~isempty(CST)
        MISSPRN = unique(setdiff(MSP(:,3),unique(CST(:,2))));
    else 
        MISSPRN = unique(MSP(:,3));
    end
    %MISSPRN'
     for prn = MISSPRN'
        MAX = [];
%         MISSDATA = sortrows(MSP,1);
        MISSDATA = sortrows(MSP(MSP(:,3)==prn,:),1);
        times = MISSDATA(:,1);
        sps = MISSDATA(:,2);
        plot(times,sps,'.-');  
        hold on;
        spmin = mean(abs(diff(sps)))+3*var(abs(diff(sps)));
%         spmin = 0.3*2*pi;
        isc = find(diff(sps)>=spmin|diff(sps)<=-spmin);
        plot(times(isc),sps(isc),'ro');
        plot(times(isc+1),sps(isc+1),'rs');
        datetick('x','HH:MM','keeplimits');
        seg = sortrows([isc;isc+1]); 
        times = times(seg);
%         [tseg,spseg] = discont_proc(times(seg),sps(seg),1/24/3600*60);
%         plot(tseg,spseg,'k.-');
%         tsegsc = (tseg-datenum([2013 12 8 3 0 0]))*24*3600;
%         tsegsc = tsegsc(~isnan(tsegsc));
        dt = 200/24/3600;
        tspan = [];
        j = 1;
        while j<=length(times)
            idx = find(times-times(j)<=dt&times-times(j)>=0);
%             if ~isempty(idx)
%                 tt = [times(j);times(idx)];
%                 tt = [tt(1);tt(end)];
%                 j = idx(end)+1;
%             else
%                 tt = [times(j);times(j)];
%                 j = j+1;
%             end
            tt = [times(j);times(idx)];
            tt = [tt(1);tt(end)];
            j = idx(end)+1;
            tspan = [tspan;tt;NaN;NaN];
        end
        odd = tspan(1:2:end);
        even = tspan(2:2:end);
        tspan_new = [odd,even];
        plot(tspan,'-');
        datetick('y','HH:MM','keeplimits');
%         datetick('x','HH:MM','keeplimits');
        %8/9
        close;
        for rr = unique(MISSDATA(:,4))'
            missdata = MISSDATA(MISSDATA(:,4)==rr,:);
            maxdata = missdata(missdata(:,2)==max(missdata(:,2)),:);
            MAX = [MAX;maxdata];
        end        
        [tt,dt] = find_tspan(MAX,span);   
        HHMM = datestr(tt,'HHMM');
%         hr_times = [hr_times;prn,str2num(year),str2num(doy),str2num(HHMM(1,:)),str2num(HHMM(2,:))];
        hr_times = [hr_times;num2str(prn,'%02i'),' ',year,' ',doy,' ',HHMM(1,:),' ',HHMM(2,:)];
        hrtimes = [hrtimes;prn tt' dt];
     end   
end
hr_times = ['PRN ' 'yr ' ' doy ' 'tstt ' 'tend';hr_times]; 
dlmwrite([op_path,'hrtimes_',year,'_',doy,'.txt'], hr_times,'delimiter','');

%read high-rate times of interest
% [hrtimes] = read_hrtimes([op_path,'hrtimes_',year,'_',doy,'.txt']);

%% High-rate processing w/{w/o} specified PRNs or time intervals
% for kk = 1:size(hrtimes(:,1),1)
prnlist = 23;
for kk = 1:size(prnlist,1)
    prn = prnlist(kk,1);
    if ismember(prn,hrtimes(:,1),'rows')
    [~,kkk] = ismember(prn,hrtimes(:,1),'rows');
    tt = hrtimes(kkk,2:3)';
    
    % specify times of interest used in the case study
    init_time = datenum([2013 12 8 3 0 0]);
    xtime = [2600;2700];
    tt = init_time + xtime/24/3600;
    disp(['Time: ',num2str(xtime(1)),'-',num2str(xtime(end)),...
        ' seconds after ',datestr(init_time,'yyyy/mm/dd HH:MM:SS')]);
%     lr_hr_sigmaphi_plot(tspan,prn,home_dir,cases_folder,rcvr_op,year,doy,sep,signal,signal_type);
    dt = diff(tt);
    
    % test with all 5 receivers including ASTRA
    disp('All 5 receivers including ASTRA');
    zoomin_hrplot_iq_carrier(tt,dt,signal_type,home_dir,cases_folder,year,doy,prn,rcvr_op);
    
    % test with 4 receivers
%     disp('4 receivers');
%     rcvr_op = ['grid108';'grid161';'grid162';'grid163'];
%     rcvr_op = 'ASTRArx';
%     zoomin_hrplot_iq_carrier(tt,dt,signal_type,home_dir,cases_folder,year,doy,prn,rcvr_op);
    end
end

end
end

