function [ tlim ] = find_tlim(datestt,dateend,doy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isempty(datestt) && ~isempty(dateend)
    d = [datestt;dateend];
    for dd = 1:2;
        xdate = d(dd,:);
        xdatestr = datestr(xdate,'mm/dd/yyyy');
        [dummy,doy_] = system(['date -d ',xdatestr,' +%j']);
        if str2num(doy_) == str2num(doy)
            ds = [xdate(1:2) xdate(3)-1 23 59 0];
            de = [xdate(1:2) xdate(3)+1 0 1 0];
        end
    end
    tlim = [datenum(ds) datenum(de)];
else
    tlim = [];

end

