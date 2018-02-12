function guicbfun(hndl)
%*************************************************************************
%*     Copyright c 2001 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
% BUGS: screws up when you open another figure because figure handle is lost
% hndl - handle of button pressed

global GUI_GPSWEEK_HNDL GUI_GPSSEC_HNDL...
       GUI_ABSSEC_HNDL...
       GUI_UTCYEAR_HNDL GUI_UTCMON_HNDL GUI_UTCDAY_HNDL GUI_UTCHOUR_HNDL GUI_UTCMIN_HNDL GUI_UTCSEC_HNDL...
       GUI_RINYEAR_HNDL GUI_RINDAY_HNDL...
       GUI_LISTBOX_HNDL GUI_WEEKDAY_HNDL;       
global GUI_COMPUTE_HNDL GUI_RESET_HNDL GUI_LEAP_HNDL;

init_gui;

if (hndl == GUI_RESET_HNDL)
    GPS.week = 0;
    GPS.sec  = 0;
    ABS.sec  = 0;
    UTC.year = 0;
    UTC.mon  = 0;
    UTC.day  = 0;
    UTC.hour = 0;
    UTC.min  = 0;
    UTC.sec  = 0;
    RINEX.year = 0;
    RINEX.day  = 0;
    set (GUI_WEEKDAY_HNDL, 'String', {'Day of week computed','0 = Sunday'});
    set (GUI_LEAP_HNDL, 'Value', 1);
elseif (hndl == GUI_COMPUTE_HNDL)
        selection = get(GUI_LISTBOX_HNDL, 'Value');
        leap = get(GUI_LEAP_HNDL, 'Value');
        switch (selection)
            case 1 % GPS Time
                GPS.week = str2num(get(GUI_GPSWEEK_HNDL, 'String'));
                GPS.sec  = str2num(get(GUI_GPSSEC_HNDL, 'String'));
                status = check_gps(GPS);
                if (status == 0)
                    [ABS, UTC, RINEX] = convert_gps(GPS, leap);
                else
                    disp('Error in GPS Time format');
                    return;
                end;
            case 2 % Absolute Time
                ABS.sec = str2num(get(GUI_ABSSEC_HNDL, 'String'));
                status = check_abs(ABS);
                if (status == 0)
                    [GPS, UTC, RINEX] = convert_abs(ABS, leap);
                else
                    disp('Error in Absolute Time format');
                    return;
                end;
            case 3 % UTC Time
                UTC.year = str2num(get(GUI_UTCYEAR_HNDL, 'String'));
                UTC.mon  = str2num(get(GUI_UTCMON_HNDL, 'String'));
                UTC.day  = str2num(get(GUI_UTCDAY_HNDL, 'String'));
                UTC.hour = str2num(get(GUI_UTCHOUR_HNDL, 'String'));
                UTC.min  = str2num(get(GUI_UTCMIN_HNDL, 'String'));
                UTC.sec  = str2num(get(GUI_UTCSEC_HNDL, 'String'));
                status = check_utc(UTC);
                if (status == 0)
                    [GPS, ABS, RINEX] = convert_utc(UTC, leap);
                else
                    disp('Error in UTC Time format');
                    return;
                end;
            case 4 % Rinex Time
                RINEX.year = str2num(get(GUI_RINYEAR_HNDL, 'String'));
                RINEX.day  = str2num(get(GUI_RINDAY_HNDL, 'String'));
                status = check_rinex(RINEX);
                if (status == 0)
                    [GPS, ABS, UTC] = convert_rinex(RINEX, leap);
                else
                    disp('Error in Rinex Time format');
                    return;
                end;
        end;
    else
    %do nothing
end;

if (hndl == GUI_RESET_HNDL | hndl == GUI_COMPUTE_HNDL)
    %Display changes
    set(GUI_GPSWEEK_HNDL, 'String', num2str(GPS.week));
    set(GUI_GPSSEC_HNDL, 'String', num2str(GPS.sec));
    
    set(GUI_ABSSEC_HNDL , 'String', num2str(ABS.sec));
    
    set(GUI_UTCYEAR_HNDL , 'String', num2str(UTC.year));
    set(GUI_UTCMON_HNDL , 'String', num2str(UTC.mon));
    set(GUI_UTCDAY_HNDL , 'String', num2str(UTC.day));
    set(GUI_UTCHOUR_HNDL , 'String', num2str(UTC.hour));
    set(GUI_UTCMIN_HNDL , 'String', num2str(UTC.min));
    set(GUI_UTCSEC_HNDL , 'String', num2str(UTC.sec));
    
    set(GUI_RINYEAR_HNDL , 'String', num2str(RINEX.year));
    set(GUI_RINDAY_HNDL , 'String', num2str(RINEX.day));
    
    days = {'Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'};
    daynum = floor(GPS.sec / (24 * 60 * 60));
    daystring = sprintf('%d = %s', daynum, char(days(daynum + 1)));
    set (GUI_WEEKDAY_HNDL, 'String', {'Day of week computed', daystring });
end;
    
    