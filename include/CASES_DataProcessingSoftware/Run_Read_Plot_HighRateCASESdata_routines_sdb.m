%% Run_Read_Plot_HighRateCASESdata_routines.m 
% % Kshitija Deshpande, Virginia Tech, May 2012. 
% Calls routines: 
% 1. Fn_ReadHighRate_CASESdata.m to read CASES high rate data
% and 
% 2. Fn_Plot_HighRate_CASESdata.m to plot CASES high rate data in several
% ways (all scintillating events, S4 and sigma_phi etc.). 
% 
% Inputs: 
% i. Needs to change "path" to enter the correct path of the CASES
% folder and folder name where all the CASES log files are stored.
% folder_path is the complete folder path including the folder name to be
% input to the read high rate data MATLAB function
% (Fn_ReadHighRate_CASESdata) and plotting high rate data function
% (Fn_Plot_HighRate_CASESdata). 
% ii. signal type:
% See Fn_ReadHighRate_CASESdata.m for details on signal type.
% Generally, it is 0 or 2 for AAL-PIP CASES. Thus, i=1 or 3. 
% iii. types of plots
% A = all & scintillating events, all segments
% of data -- processed
% B = S4 and Sigma_phi plots
% C = plot common scintillating times for PRNs
%
% Adapting to use at IIT. S. Datta-Barua, 11 Sept 2013.


close all
% clear all

sep = filesep;

%period in question;
tstart = [2013 6 23 10 34 0];
tend = [2013 6 23 10 39 0];
hour = tstart(4);

%User defined path to the log files of the data.
%Change it for your system
%path = '/Users/kshitija/work/research/CASES_Data/2012/'; 
% path = '/Users/kshitija/work/research/CASES_Data/10-dayCampaignNov2012/'; 
case_folder = '/home/yang/cases/';

%Change folder name according to the date of CASES data, folder with all
%the log files
% folder_name = 'ComparingWith2012_01_24/2012_01_24';
rcvr_name = 'grid108';
doy = '174';
folder_path = [case_folder,rcvr_name,sep,doy,sep];

disp(folder_path)

for i = 1:1:1; % for signal type,
    %see Fn_ReadHighRate_CASESdata.m for details on signal type
    signal_type = i-1;
    Fn_ReadHighRate_CASESdata_sdb(folder_path,signal_type,doy,tstart,tend,hour);
    for ii = 1:1:1; %Type of plots
        %specify types of plots
        % ii = 1 -> A = all & scintillating events, all segments
        % of data -- processed
        % ii = 2 -> B = S4 and Sigma_phi plots
        % ii = 3 -> C = plot common scintillating times for PRNs
        set_plotCT = ii-1;
        switch set_plotCT
            case 0
                set_plot = 'A';
            case 1
                set_plot = 'B';
            case 2
                set_plot = 'C';
            otherwise
                error('Unknown set.')
        end
        Fn_Plot_HighRate_CASESdata_sdb(folder_path,signal_type,set_plot);
    end %set_plot
end %signal type
