function status = check_gps(GPS)

% Check the validity of Data in the struct GPS
% GPS.week
% GPS.sec
% status of 0 implies data is fine, anything else is an error

status = 0;

if (isempty(GPS.week) == 1 | isempty(GPS.sec) == 1) %wasn't a number
    status = 1;
    disp('GPS Time information not numeric');
    return;
end;

if (GPS.week - ceil(GPS.week) ~= 0.0)
    status = 1;
    disp('GPS week must be a whole number');
    return;
end;

if (GPS.sec >= 60*60*24*7)
    status = 1;
    disp('GPS seconds cannot be greater than 1 week long');
    return;
end;

if (GPS.sec < 0 | GPS.week < 0)
    status = 1;
    disp('GPS Time information cannot be negative');
end;

return;