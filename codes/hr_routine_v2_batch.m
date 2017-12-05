function [] = hr_routine_v2_batch(yearin, doyin)
%% Initialization
close all;
dbstop if error;
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
[~, yesterday] = system(comm);
yesterday = str2num(yesterday);

comm1 = 'date -u -d "a day ago" +%Y';
[~, year] = system(comm1);

if nargin == 2
    yearlist = yearin;
    doylist = doyin;
else
    yearlist = str2num(year);
    doylist = yesterday;
end

SD = [];
MEGA_MSP = [];
MSP_days = [];
SCINTEVENTS = [];
%specify signal
for signal_type = [0, 2]
    for yearnum = yearlist
        year = num2str(yearnum, '%04i')
        % for doy = yesterday
        for doynum = unique(doylist)
            %day of year in string
            doy = num2str(doynum, '%03i')
            year;
            %% Process low-rate data
            %first sigmaphi threshold in [rad]
            spmask = 0; s4mask = 0;
            
            if strcmp(cases_folder(end-4:end-1), 'pfrr')
                %folder_path for 2013 Poker Flat data
                op_path = [home_dir, 'PFRR_Data/'];
            else
                %folder_path for 2013 Calgary data
                op_path = [home_dir, 'Calgary_Data/'];
            end
            
            %check if low rate data has already been processed
            matfile = dir([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'])
            if ~isempty(matfile) && matfile.datenum >= datenum([2015, 5, 30, 0, 0, 0])
                disp([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat already exists']);
            else
                %filter low rate sigma phi data for all operational receivers
                %with sigmaphi threshold and save the results into a large array
                [MSP, MS4, HITDATA, CST, rcvr_op, tlim, splim, s4lim, signal] = scint_el_stackplot(cases_folder, home_dir, signal_type, doy, year, spmask, s4mask);
                lrdata = {MSP, HITDATA, CST, rcvr_op, tlim, splim, signal, spmask, MS4, s4mask, s4lim};
                if ~isempty(MSP) && ~isempty(MS4)
                    save([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'], 'lrdata');
                else
                    disp('No data for this day, exiting....');
                    continue;
                end
            end
            %% Generate low-rate sigmaphi stackplots
            %load low rate results
            lrfile = [op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'];
            load(lrfile);
            tlim_vec = datevec(lrdata{5})
            hr_times = [];
            hrtimes = [];
            MSP = lrdata{1}; HITDATA = lrdata{2};
            CST = lrdata{3}; rcvr_op = lrdata{4};
            tlim = lrdata{5}; splim = lrdata{6};
            signal = lrdata{7}; spmask = lrdata{8};
            MS4 = lrdata{9}; s4mask = lrdata{10};
            s4lim = lrdata{11};
            span = 60;
            
            disp(['There are totally ', num2str(size(MSP, 1)), ...
                ' points of sigma_phi data']);
            disp(['There are totally ', num2str(size(MS4, 1)), ...
                ' points of s4 data']);
            
            % receiver structure including all poker flat receivers
            rcvr_struct = rx_dirs(cases_folder, year, doy);
            rcvr_struct = ['grid108'; 'grid154'; 'grid160'; 'grid161'; 'grid162'; 'grid163'; 'ASTRArx'];
            % rcvr_struct = ['grid108';'grid154';'grid160';'grid161';'grid162';'grid163'];
            
            if ~isempty(MSP)
                lr_rx_stackplot(lrdata, year, doy, rcvr_struct, op_path, 'sp');
                lr_prn_stackplot(lrdata, year, doy, rcvr_struct, op_path, 'sp');
            else % no enough scintillation data to conitue, considered a quiet day
                disp(['doy:', doy, ' of ', year, ' is a quiet day']);
            end
            
            if ~isempty(MS4)
                lr_rx_stackplot(lrdata, year, doy, rcvr_struct, op_path, 's4');
                lr_prn_stackplot(lrdata, year, doy, rcvr_struct, op_path, 's4');
            else % no enough scintillation data to conitue, considered a quiet day
                disp(['doy:', doy, ' of ', year, ' is a quiet day']);
                continue;
            end
        end
    end
end
end