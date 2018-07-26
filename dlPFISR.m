clear all;
% globalDownload('http://isr.sri.com/madrigal/', ...
%     '.', ...
%     'Yang Su', ...
%     'ysu27@hawk.iit,edu', ...
%     'IIT', ...
%     'hdf5', ...
%     datenum('17-Mar-2015 00:00:00'), ...
%     datenum('18-Mar-2015 23:59:59'), ...
%     61, ...
%     5960, ...
%     '', ...
%     '')

globalIsprint('http://isr.sri.com/madrigal', ...
    'YEAR,MONTH,DAY,HOUR,MIN,SEC,INTTMS,CGM_LAT,VIPE,VIPN,VI6,DVIPE,DVIPN,DVI6,MAGVEL,NANGLE,DMAGVEL,DNANGLE', ...
    'PFISR_data1.txt', ...
    'Yang Su', ...
    'ysu27@hawk.iit.edu', ...
    'IIT', ...
    datenum('03-17-2015 00:00:00'), ...
    datenum('03-19-2015 23:59:59'), ...
    61, ... %instrument code for PFISR
    'filter=CGM_LAT,,70', ...
    5960, ... %vector velocity
    '', ...
    '')

% globalIsprint('http://isr.sri.com/madrigal', ...
%     'YEAR,MONTH,DAY,HOUR,MIN,SEC,AZM,ELM,BEAMID,RANGE,NEL,DNEL,TI,DTI,TE,DTE', ...
%     'PFISR_data2.txt', ...
%     'Yang Su', ...
%     'ysu27@hawk.iit.edu', ...
%     'IIT', ...
%     datenum('12-08-2013 00:00:00'), ...
%     datenum('12-08-2013 23:59:59'), ...
%     61, ... %instrument code for PFISR
%     'filter=CGM_LAT,,70', ...
%     5950, ... %long pulse
%     '', ...
%     '')