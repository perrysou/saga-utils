function el_az_stackplot(folder_path,signal_type,rcvr_name)
close all;
clear all;
folder_path = 'D:\MMAE552\~sdattaba\downloads\grid108\173\';
% folder_path = 'C:\cases\grid108\173\';
signal_type = 0;
rcvr_name = 'grid108';
sep = filesep;
%specify signal type
switch signal_type
   case 0
      signal = 'L1CA';  %// GPS L1 legacy civil code 
   case 1
      signal = 'L2CM';  %// GPS L2 civil L code
   case 2
      signal = 'L2CL';  %// GPS L2 M+L combined tracking
   case 3
      signal = 'L2CLM';  %// GPS L5 civil in-phase
   case 4
      signal = 'L5I';   %// GPS L5 civil quadrature
   case 5
      signal = 'L5Q';
   case 6
      signal = 'L5IQ';
   case 7
      signal = 'L1CA-ALT1';
   case 8
      signal = 'CDMA-UHF-PILOT';
   case 9
      signal = 'CDMA-UHF-SYNC';
   otherwise
      error('Unknown signal.')
end

%Create output folders
command = strcat('mkdir',{' '},folder_path,'plots_',signal);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_path,'plots_',signal,sep,'LowRate',sep);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_path,'prn_files_',signal);
system(cell2mat(command));
outdir = strcat(folder_path,'plots_',signal,sep,'LowRate',sep);
sep = filesep;

%Read and combine txfinto files
azelfilestruct = strcat(folder_path,'txt',sep,ls([folder_path,'txt',sep,'txinfo*.log']));
AZEL=[];
for jj=1:length(azelfilestruct(:,1))
    azelfile = dlmread(azelfilestruct(jj,:));  
    AZEL = [AZEL;azelfile];
end

%Remove invalid entries data
AZEL = AZEL((AZEL(:,5)>0),:);
AZEL_old = AZEL;
%Low elevation mask;
AZEL_low = AZEL(AZEL(:,5)>30,:);
% AZEL_low = sortrows(AZEL_low,8);

%Gather data
ORTW = AZEL(:,1);
ORTS = AZEL(:,2)+AZEL(:,3);
EL = AZEL(:,5);
PRN = AZEL(:,8);
GPSTIME =  [ORTW ORTS];
L_GPSTIME = size(GPSTIME,1);
UTC_init = gps2utc(AZEL(1,1:2),0);
UTC_end = gps2utc(AZEL(end,1:2),0);
YEAR_init = UTC_init(1);
MONTH_init = UTC_init(2);
DAY_init = UTC_init(3);
HOUR_init = UTC_init(4);
HOUR_end = UTC_end(4);

%%
DATA_el = [ORTW ORTS EL PRN];
% DATA = sortrows(DATA,4);
outfilename = strcat('prn_files_',signal,sep,'txinfodata.mat');
save([folder_path,outfilename],'DATA_el');
DATA = DATA_el;

%Create mega info arrays for every PRN
KK = unique(DATA(:,4));
for kk = KK(1):KK(end)
    sv = find(DATA(:,4)==kk);
    svdata = DATA(sv,:);
    ortw = DATA(sv,1);
    orts = DATA(sv,2); 
    obs_time = gps2utc([ortw orts],0);
    el = DATA(sv,3);
    if size(obs_time,1)>4
        time = (obs_time(:,4)-HOUR_init) + obs_time(:,5)/60 + obs_time(:,6)/3600;
        dt = min(diff(time));
        t_d = find(diff(time)>0.5)';
        seq_d = [1 t_d length(time)];
        time_e = [];el_e = [];
        if isempty(t_d)==1
            time_e = time;
            el_e = el;
        else
            for ind_d = 1:1:length(t_d)
                time_d = (time(seq_d(ind_d+1))+dt:dt:time(seq_d(ind_d+1)+1))';
                time_d = time_d*NaN';
                el_d = ones(size(time_d))*NaN;
                if ind_d~=1
                   time_e = [time_e;time(seq_d(ind_d)+1:seq_d(ind_d+1));time_d];
                    el_e = [el_e;el(seq_d(ind_d)+1:seq_d(ind_d+1));el_d];
                else
                    time_e = [time_e;time(seq_d(ind_d):seq_d(ind_d+1));time_d];
                    el_e = [el_e;el(seq_d(ind_d):seq_d(ind_d+1));el_d];
                end
            end
            time_e = [time_e;time(seq_d(end-1)+1:seq_d(end))];
            el_e = [el_e;el(seq_d(end-1)+1:seq_d(end))];
        end
    subplot(8,4,kk);
    plot(time_e,el_e,'LineWidth',1.5,'Color',rand(1,3));
    xlim([HOUR_init HOUR_end+1]);ylim([0 90]);
    xlabel(['PRN',num2str(kk)]);
    set(gca,'XTick',HOUR_init:2:HOUR_end+1,'YTick',15:15:30,'YGrid','on','XGrid','off');
    end
end
end
    


