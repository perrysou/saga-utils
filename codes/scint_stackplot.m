function scint_stackplot(folder_path,signal_type,rcvr_name)
close all;
clear all;
folder_path = 'D:\MMAE552\~sdattaba\downloads\grid160\173\';
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
scintfilestruct = strcat(folder_path,'txt',sep,ls([folder_path,'txt',sep,'scint*.log']));
SCINT=[];
for jj=1:length(scintfilestruct(:,1))
    scintfile = dlmread(scintfilestruct(jj,:));  
    SCINT = [SCINT;scintfile];
end

%Specify signal type
SCINT = SCINT((SCINT(:,14)==signal_type),:);
%Remove invalid entries
SCINT = SCINT((SCINT(:,8)<5)&SCINT(:,8)>-1,:);
SCINTNS = unique(SCINT((SCINT(:,13)==1),15))

%Gather data
ORTW = SCINT(:,1);
ORTS = SCINT(:,2)+SCINT(:,3);
SIGMAPHI = SCINT(:,8);
PRN = SCINT(:,15);
GPSTIME =  [ORTW ORTS];
L_GPSTIME = size(GPSTIME,1);
UTC_init = gps2utc(SCINT(1,1:2),0);
UTC_end = gps2utc(SCINT(end,1:2),0);
YEAR_init = UTC_init(1);
MONTH_init = UTC_init(2);
DAY_init = UTC_init(3);
HOUR_init = UTC_init(4);
HOUR_end = UTC_end(4);
UTC_TIME = gps2utc(GPSTIME,0);
TIME = (UTC_TIME(:,4)-HOUR_init) + UTC_TIME(:,5)/60 + UTC_TIME(:,6)/3600;

%%
DATA_scint = [GPSTIME SIGMAPHI PRN];
% DATA = sortrows(DATA,4);
outfilename = strcat('prn_files_',signal,sep,'scintdata.mat');
save([folder_path,outfilename],'DATA_scint');

DATA = DATA_scint;
%Create mega info arrays for every PRN
KK = unique(DATA(:,4));

load D:\MMAE552\~sdattaba\downloads\grid108\173\prn_files_L1CA\txinfodata.mat

SC = [];
for kk = KK(1):KK(end)
    sv = find(DATA(:,4)==kk);
    svdata = DATA(sv,:);
    ortw = DATA(sv,1);
    orts = DATA(sv,2); 
    time_sc = gps2utc([ortw orts],0);
    time_sc = (time_sc(:,4)-HOUR_init) + time_sc(:,5)/60 + time_sc(:,6)/3600;
    sv_el = find(DATA_el(:,4)==kk & DATA_el(:,3)>15);
    ortw_el = DATA_el(sv_el,1);  
    orts_el = DATA_el(sv_el,2);   
    if isempty(orts_el) || isempty(ortw_el)== 1
        disp(['PRN ',num2str(kk),' is at low elevation(<15 degrees).']);
    else
        time_el = gps2utc([ortw_el orts_el],0);
        time_el = (time_el(:,4)-HOUR_init) + time_el(:,5)/60 + time_el(:,6)/3600;
        svdata_masked = svdata(time_sc<=max(time_el) & time_sc>=min(time_el),:);
        ortw_masked = svdata_masked(:,1);orts_masked = svdata_masked(:,2);
        sigmaphi = svdata_masked(:,3);
        obs_time = gps2utc([ortw_masked orts_masked],0);
        if size(obs_time,1)>4
            time = (obs_time(:,4)-HOUR_init) + obs_time(:,5)/60 + obs_time(:,6)/3600;
            dt = min(diff(time));
            t_d = find(diff(time)>0.5)';
            seq_d = [1 t_d length(time)];
            time_e = [];sigmaphi_e = [];
            if isempty(t_d)==1
                time_e = time;
                sigmaphi_e = sigmaphi;
            else
                for ind_d = 1:1:length(t_d)
                    time_d = (time(seq_d(ind_d+1))+dt:dt:time(seq_d(ind_d+1)+1))';
                    time_d = time_d*NaN';
                    sigmaphi_d = ones(size(time_d))*NaN;
                    if ind_d~=1
                       time_e = [time_e;time(seq_d(ind_d)+1:seq_d(ind_d+1));time_d];
                        sigmaphi_e = [sigmaphi_e;sigmaphi(seq_d(ind_d)+1:seq_d(ind_d+1));sigmaphi_d];
                    else
                        time_e = [time_e;time(seq_d(ind_d):seq_d(ind_d+1));time_d];
                        sigmaphi_e = [sigmaphi_e;sigmaphi(seq_d(ind_d):seq_d(ind_d+1));sigmaphi_d];
                    end
                end
                time_e = [time_e;time(seq_d(end-1)+1:seq_d(end))];
                sigmaphi_e = [sigmaphi_e;sigmaphi(seq_d(end-1)+1:seq_d(end))];
            end
            subplot(8,4,kk);
            set(gca,'FontSize',6);
            plot(time_e,sigmaphi_e,'LineWidth',1.5,'Color',rand(1,3));
            xlim([HOUR_init HOUR_end+1]);ylim('auto');
            title(strcat('Low Rate \sigma_\phi for PRN',num2str(kk)));
            grid on;
            set(gca,'XTick',HOUR_init:2:HOUR_end+1);
            sctime = time_e(sigmaphi_e==max(sigmaphi_e),:);
            SC = [SC; sctime max(sigmaphi_e)*ones(size(sctime)) kk*ones(size(sctime))]
        end
    end
    
end
% figure = gca;
% title = strcat({'Low rate standard deviation of phase, \sigma_{\phi}'});
fig_name = strcat(signal,'_LowRateSigmaPhi');
fig = strcat(outdir,fig_name,'.fig');
saveas(gcf,fig);
    


