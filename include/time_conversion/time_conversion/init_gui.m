% globals for GUI objects

global GUI_GPSWEEK_HNDL GUI_GPSSEC_HNDL...
       GUI_ABSSEC_HNDL...
       GUI_UTCYEAR_HNDL GUI_UTCMON_HNDL GUI_UTCDAY_HNDL GUI_UTCHOUR_HNDL GUI_UTCMIN_HNDL GUI_UTCSEC_HNDL...
       GUI_RINYEAR_HNDL GUI_RINDAY_HNDL...
       GUI_LISTBOX_HNDL...
       GUI_WEEKDAY_HNDL;
global GUI_COMPUTE_HNDL GUI_RESET_HNDL GUI_LEAP_HNDL;

listbox = {'GPS Time', 'Absolute Time', 'UTC Time', 'Rinex File Time'};

GUI_GPSWEEK_HNDL = findobj('Tag', 'GPS_week');
GUI_GPSSEC_HNDL =  findobj('Tag', 'GPS_sec');

GUI_ABSSEC_HNDL =  findobj('Tag', 'ABS_sec');

GUI_UTCYEAR_HNDL = findobj('Tag', 'UTC_year');
GUI_UTCMON_HNDL =  findobj('Tag', 'UTC_mon');
GUI_UTCDAY_HNDL =  findobj('Tag', 'UTC_day');
GUI_UTCHOUR_HNDL = findobj('Tag', 'UTC_hour');
GUI_UTCMIN_HNDL =  findobj('Tag', 'UTC_min');
GUI_UTCSEC_HNDL =  findobj('Tag', 'UTC_sec');

GUI_RINYEAR_HNDL = findobj('Tag', 'RIN_year');
GUI_RINDAY_HNDL =  findobj('Tag', 'RIN_day');

GUI_LISTBOX_HNDL = findobj('Tag', 'input_list');
set(GUI_LISTBOX_HNDL, 'String', listbox);
GUI_RESET_HNDL = findobj('Tag','Reset');
GUI_COMPUTE_HNDL = findobj('Tag','Calculate');
GUI_WEEKDAY_HNDL = findobj('Tag', 'weekday');
GUI_LEAP_HNDL = findobj('Tag', 'leap_seconds');














