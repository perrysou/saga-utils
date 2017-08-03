function [] = hr_routine_v1_batch(yearin,doyin)
% hr_rooutine_v1_batch
% used only for batch processing
%hr_routine_v1_batch(2013,1:366)
% 3/27/2015 Yang Su
%% Initialization
close all;
warning off;
% dbstop if error;
format long g;
% maxNumCompThreads(4);

%file separator "/" in linux
sep = filesep;

%where cases data are stored
cases_folder = '/data1/public/Data/cases/pfrr/';
%case_folder = '/data1/public/Data/cases/calg/';

%where output data and plots are stored
home_dir = '/data1/home/ysu27/';

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

%specify signal 
for signal_type = 0 %[0 2]
% for doy = yesterday
for doynum = doylist
%day of year in string
doy = num2str(doynum,'%03i')
year
%% Process low-rate data 
%first sigmaphi threshold in [rad]
spmask = 0;s4mask = 0;
%filter low rate sigma phi data for all operational receivers
%with sigmaphi threshold and save the results into a large array
[MSP,MS4,HITDATA,CST,rcvr_op,tlim,splim,s4lim,signal] = scint_el_stackplot(cases_folder,home_dir,signal_type,doy,year,spmask,s4mask);
lrdata = {MSP,HITDATA,CST,rcvr_op,tlim,splim,signal,spmask,MS4,s4mask,s4lim};
if strcmp(cases_folder(end-4:end-1),'pfrr')
    %folder_path for 2013 Poker Flat data 
    op_path = [home_dir,'PFRR_Data/'];
else
    %folder_path for 2013 Calgary data
    op_path = [home_dir,'Calgary_Data/'];  
end
save([op_path,'lrdata_',year,'_',doy,'.mat'],'lrdata');
%% Generate low-rate sigmaphi stackplots 
%load low rate results
load([op_path,'lrdata_',year,'_',doy,'.mat']);
% load('home/yang/Dropbox/research/lrdata_2013_342.mat');
hr_times = [];
hrtimes = [];
MSP = lrdata{1};HITDATA = lrdata{2};
CST = lrdata{3};rcvr_op = lrdata{4};
tlim = lrdata{5};splim = lrdata{6};
signal = lrdata{7};spmask = lrdata{8};
MS4 = lrdata{9};s4mask = lrdata{10};
s4lim = lrdata{11};
span = 60;

disp(['There are totally ', num2str(size(MSP,1)), ...
    ' points of sigma_phi data after the first filter']);
disp(['There are totally ', num2str(size(MS4,1)), ...
    ' points of s4 data after the first filter']);

% receiver structure including all poker flat receivers
rcvr_struct = rx_dirs(cases_folder,year,doy);
rcvr_struct = ['grid108';'grid154';'grid160';'grid161';'grid162';'grid163';'ASTRArx'];

if ~isempty(MSP)
    lr_rx_stackplot(lrdata,year,doy,rcvr_struct,op_path,'sp');
    lr_prn_stackplot(lrdata,year,doy,rcvr_struct,op_path,'sp');
else % no enough scintillation data to conitue, considered a quiet day
    disp(['doy:',doy,' of ',year,' is a quiet day']);
end

if ~isempty(MS4)
    lr_rx_stackplot(lrdata,year,doy,rcvr_struct,op_path,'s4');
    lr_prn_stackplot(lrdata,year,doy,rcvr_struct,op_path,'s4');
else % no enough scintillation data to conitue, considered a quiet day
    disp(['doy:',doy,' of ',year,' is a quiet day']);
end

continue;
end
end
end

