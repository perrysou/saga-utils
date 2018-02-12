function status = check_rinex(RINEX)

% Check the validity of Data in the struct RINEX
% RINEX.year
% RINEX.day
% status of 0 implies data is fine, anything else is an error

status = 0;

if (isempty(RINEX.year) | RINEX.year < 1980)
    status = 1;
    return;
    disp('Invalid Rinex file year');
end;

if (isempty(RINEX.day) | RINEX.day < 1 | RINEX.day > 366)
    status = 1;
    return;
    disp('Invalid RINEX file day');
end;

if (RINEX.year - ceil(RINEX.year) ~= 0.0 |...
    RINEX.day - ceil(RINEX.day) ~= 0.0)
    status = 1;
    disp('Rinex file data not a whole number');
    return;
end;

if (isleapyear(RINEX.year) == 0 & RINEX.day == 366 |...
   (RINEX.year == 1980 & RINEX.day < 6))
    status = 1;
    disp('Invalid Rinex File year-day combination');
        return;
end;

return;