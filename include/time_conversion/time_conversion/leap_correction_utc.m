function delta = leap_correction_utc(UTC)

% return number of leap seconds included at time UTC
% where UTC is a structure containing year, mon, day, hour, min,
% sec.

if (UTC.year < 1981 | (UTC.year == 1981 & UTC.mon < 7)) % July 1, 1981
    delta = 0;
elseif(UTC.year < 1982 | (UTC.year == 1982 & UTC.mon < 7)) % July 1, 1982
    delta = 1;
elseif(UTC.year < 1983 | (UTC.year == 1983 & UTC.mon < 7)) % July 1, 1983
    delta = 2;
elseif(UTC.year < 1985 | (UTC.year == 1985 & UTC.mon < 7)) % July 1, 1985
    delta = 3;
elseif(UTC.year < 1988) % Jan 1, 1988
    delta = 4;
elseif(UTC.year < 1990) % Jan 1, 1990
    delta = 5;
elseif(UTC.year < 1991) % Jan 1, 1991
    delta = 6;
elseif(UTC.year < 1992 | (UTC.year == 1992 & UTC.mon < 7)) % July 1, 1992
    delta = 7;
elseif(UTC.year < 1993 | (UTC.year == 1993 & UTC.mon < 7)) % July 1, 1993
    delta = 8;
elseif(UTC.year < 1994 | (UTC.year == 1994 & UTC.mon < 7)) % July 1, 1994
    delta = 9;
elseif(UTC.year < 1996) % Jan 1, 1996
    delta = 10;
elseif(UTC.year < 1997 | (UTC.year == 1997 & UTC.mon < 7)) % July 1, 1997
    delta = 11;
elseif(UTC.year < 1999) % Jan 1, 1999
    delta = 12;
else
    delta = 13;
end;