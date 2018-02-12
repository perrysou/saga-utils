function status = check_utc(UTC)

% Check the validity of Data in the struct UTC
% UTC.year
% UTC.mon
% UTC.day
% UTC.hour
% UTC.min
% UTC.sec
% status of 0 implies data is fine, anything else is an error

status = 0;

if (isempty(UTC.year) | UTC.year < 1980)
    status = 1;
    disp('Invalid UTC year');
    return;
end;

if (isempty(UTC.mon) | UTC.mon < 1 | UTC.mon > 12)
    status = 1;
    disp('Invalid UTC month');
    return;
end;

if (isempty(UTC.day) | UTC.day < 1 | UTC.day > 31)
    status = 1;
    disp('Invalid UTC day');
    return;
end;

if (isempty(UTC.hour) | UTC.hour < 0 | UTC.hour > 23)
    status = 1;
    disp('Invalid UTC hour');
    return;
end;

if (isempty(UTC.min) | UTC.min < 0 | UTC.min > 59)
    status = 1;
    disp('Invalid UTC minute');
    return;
end;

if (isempty(UTC.sec) | UTC.sec < 0 | UTC.sec >= 60)
    status = 1;
    disp('Invalid UTC second');
    return;
end;

if (UTC.year - ceil(UTC.year) ~= 0.0 |...
    UTC.mon - ceil(UTC.mon) ~= 0.0 |...
    UTC.day - ceil(UTC.day) ~= 0.0 |...
    UTC.hour - ceil(UTC.hour) ~= 0.0 |...
    UTC.min - ceil(UTC.min) ~= 0.0)
    status = 1;
    disp('UTC data not a whole number');
    return;
    return;
end;

if ((UTC.mon == 4 & UTC.day == 31) |...
    (UTC.mon == 6 & UTC.day == 31) |...
    (UTC.mon == 9 & UTC.day == 31) |...
    (UTC.mon == 11 & UTC.day == 31) |...
    (UTC.mon == 2 & UTC.day > 29) |...
    (UTC.mon == 2 & isleapyear(UTC.year) == 0 & UTC.day == 29) |...
    (UTC.year == 1980 & UTC.mon == 1 & UTC.day < 6))
        status = 1;
        disp('Invalid Year-Month-Day combination');
        return;
end;

return;