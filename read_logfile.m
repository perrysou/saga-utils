function [file] = read_logfile(filetype, tstart, tend, rcvr_name, cases_folder, op_path)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%---------
%test
% filetype = 'iq*';
% tstart = [2013 12 8 22 43 20];
% tend = [2013 12 9 2 45 0];
% rcvr_name = 'grid108';
% cases_folder = '/data1/public/Data/cases/pfrr/';
% op_path = '/data1/home/ysu27/PFRR_Data/grid108/2013/342/';
%test end
%---------

global sep
%remove xdata.log from previous aborted data processing
% comm = ['find ',op_path,' -name ''xdata.log'' -exec rm {} \;'];
comm = ['rm -v ', op_path, 'xdata.log'];
system(comm);

hourstt = tstart(4);
hourend = tend(4);
yearstt = tstart(1);
yearend = tend(1);
%find days of year for start time and end time
tstt_str = datestr(tstart, 'mm/dd/yyyy');
tend_str = datestr(tend, 'mm/dd/yyyy');
if strcmp('MACI64', computer)
    [~, doystt] = system(['date -ujf "%m/%d/%Y" ', tstt_str, ' +%j']);
    [~, doyend] = system(['date -ujf "%m/%d/%Y" ', tend_str, ' +%j']);
end
if strcmp('GLNXA64', computer)
    [~, doystt] = system(['date -ud ', tstt_str, ' +%j']);
    [~, doyend] = system(['date -ud ', tend_str, ' +%j']);
end
doystt = str2num(doystt);
doyend = str2num(doyend);
% find gps times(gps week and second) for start time and end time
tspan_gps = utc2gps([tstart; tend], 0);
gpsstt = tspan_gps(1, :);
gpsend = tspan_gps(2, :);
gpswstt = gpsstt(1);
gpssstt = gpsstt(2);
gpswend = gpsend(1);
gpssend = gpsend(2);

cases_year = [cases_folder, num2str(yearstt, '%04i'), sep];
DATA = [];
if strcmp(filetype, 'iq*') == 1 %read iqdata
    %find offsetdata
    if hourstt == 0
        offsett_str = datestr(datenum(tstart)-3600/24/3600, 'mm/dd/yyyy HH:MM:SS');
        [~, offsett] = system(['date -ud ''', offsett_str, ''' +%Y%j%H']);
        offsetyear = str2num(offsett(1:4));
        offsetdoy = str2num(offsett(5:7));
        offsethour = str2num(offsett(end-2:end-1));
        sub_read_logfile(offsetyear, offsetdoy, sep, rcvr_name, offsethour, filetype, ...
            cases_folder, cases_year, op_path, gpswstt, gpssend, gpssstt);
        hourstt = hourstt + 1;
        % add lines to read offset data;
    end
    if hourstt > hourend
        hour = hourstt - 1;
        while hour <= 23
            hour = hour + 1;
        end
        hour = hour - 24;
        while hour <= hourend + 1
            hour = hour + 1;
        end
    elseif hourstt <= hourend %mostly used condition, refer here for other conditions
        hour = hourstt - 1;
        while hour <= hourend + 1
            sub_read_logfile(yearstt, doyend, sep, rcvr_name, hour, filetype, ...
                cases_folder, cases_year, op_path, gpswstt, gpssend, gpssstt);
            hour = hour + 1;
        end
    end
    
    
end
% !head -n 1 /data1/home/ysu27/PFRR_Data/grid108/2014/019/xdata.log
% !tail -n 1 /data1/home/ysu27/PFRR_Data/grid108/2014/019/xdata.log
if ~isempty(dir([op_path, 'xdata.log']))
    file = importdata([op_path, 'xdata.log']);
else
    disp({'Caution! You are not supposed to get to this step as operational receivers are determined beforehand'; ...
        'this receiver may be unoperational'});
    file = [];
    comm = ['rm -v ', op_path, 'xdata.log'];
    system(comm);
    % %repack iq logfiles
    % repack_comm = ['sh ',pwd,sep,'repack_v1.sh ',in_path_tmp,' ',logfile_path,' ',filetype]
    % system(repack_comm);
end
end

function [] = sub_read_logfile(yearstt, doyend, sep, rcvr_name, hour, filetype, ...
    cases_folder, cases_year, op_path, gpswstt, gpssend, gpssstt)
global home_dir;
%iq binaries in /data2/from_usb/
in_path0 = strcat(['/data2/from_usb/', num2str(yearstt, '%04i'), sep, num2str(doyend, '%03i'), sep, rcvr_name, sep]);
flag0 = dir([in_path0, 'bin', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.bin']);

%iq binaries downloaded on demand in /data1/public/Data/cases/pfrr/from_usb
in_path1 = strcat([cases_folder, sep, 'from_usb', sep, num2str(yearstt, '%04i'), sep, num2str(doyend, '%03i'), sep, rcvr_name, sep]);
flag1 = dir([in_path1, 'bin', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.bin']);

%iq binaries in /data2/from_usb2/
in_path3 = strcat(['/data2/from_usb2/', num2str(yearstt, '%04i'), sep, num2str(doyend, '%03i'), sep, rcvr_name, sep]);
flag3 = dir([in_path3, 'bin', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.bin']);

%iq binaries in /data2/from_usb3/
in_path4 = strcat(['/data2/from_usb3/', num2str(yearstt, '%04i'), sep, num2str(doyend, '%03i'), sep, rcvr_name, sep]);
flag4 = dir([in_path4, 'bin', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.bin']);

%iq binaries already downloaded in /data1/public/Data/cases/pfrr/
in_path2 = strcat([cases_year, num2str(doyend, '%03i'), sep, rcvr_name, sep]);
flag2 = dir([in_path2, 'bin', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.bin']);
[in_, ~] = inoutpath(cases_folder, home_dir, ...
    num2str(yearstt, '%04i'), num2str(doyend, '%03i'), rcvr_name);

if strcmp(rcvr_name, 'ASTRArx')
    flag2 = dir([in_, 'txt', sep, '*', filetype, '_', num2str(hour, '%02i'), '*.log']);
    %     keyboard;
end
if isempty(flag0) && isempty(flag1) && isempty(flag2) && isempty(flag3) && isempty(flag4)
    disp(['Unable to find any binaries for hr:', num2str(hour, '%02i'), ' on the server, please download manually']);
    %exception for IIT-13 data for doy 342, 2013 since only logfiles are
    %available
    if ~strcmp(rcvr_name, 'ASTRArx')
        return;
    end
elseif ~isempty(flag2)
    in_path = in_path2;
elseif ~isempty(flag0)
    in_path = in_path0;
elseif ~isempty(flag1)
    in_path = in_path1;
elseif ~isempty(flag3)
    in_path = in_path3;
elseif ~isempty(flag4)
    in_path = in_path4;
end

logfile_dir = [op_path(1:end-9), num2str(yearstt, '%04i'), sep, num2str(doyend, '%03i'), sep];
logfile_name = [logfile_dir, 'txt', sep, filetype, '_', num2str(yearstt, '%04i'), ...
    '_', num2str(doyend, '%03i'), '_', num2str(hour, '%02i'), '*.log'];
if isempty(dir(logfile_name))
    if strcmp(rcvr_name, 'ASTRArx')
        %         return;
    else
        unpack_comm = ['sh ', '/data1/home/ysu27', sep, 'unpack_v1.sh ', in_path, ' ', logfile_dir, ' ', filetype, num2str(hour, '%02i'), '??'];
        system(unpack_comm);
    end
else
    disp(['iq binaries already unpacked; locate logfiles in ', logfile_dir, 'txt', sep]);
end
if ~strcmp(rcvr_name, 'ASTRArx')
    comm = ['awk ''$3==', num2str(gpswstt), ...
        '&&$4<=', num2str(gpssend), '&&$4>=', num2str(gpssstt), ...
        ' {print $0;}'' ', logfile_dir, 'txt', sep, ...
        filetype, num2str(hour, '%02i'), '??.log', ' >> ', op_path, 'xdata.log'];
else
    comm = ['awk ''$3==', num2str(gpswstt), ...
        '&&$4<=', num2str(gpssend), '&&$4>=', num2str(gpssstt), ...
        ' {print $0;}'' ', in_, 'txt', sep, ...
        filetype, num2str(hour, '%02i'), '??.log', ' >> ', op_path, 'xdata.log'];
end
system(comm);
end