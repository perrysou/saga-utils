function status = isleapyear(year)

% return 1 if this year is a leap year
% return 0 if it isn't

status = 0;

if (mod(year, 4) == 0 )
    if ( mod(year, 100) == 0 & mod(year, 400) > 0) % Little known rule for leap years
        return;
    end;
    status = 1;
end;

return;