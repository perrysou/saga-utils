function delta = leap_correction_gps(total)

% return number of leap seconds not included at time 'total'
% where total is the number of seconds since Jan 6, 1980

if (total < 46828800) % July 1, 1981
    delta = 0;
elseif(total < 78364800) % July 1, 1982
    delta = 1;
elseif(total < 109900800) % July 1, 1983
    delta = 2;
elseif(total < 173059200) % July 1, 1985
    delta = 3;
elseif(total < 252028800) % Jan 1, 1988
    delta = 4;
elseif(total < 315187200) % Jan 1, 1990
    delta = 5;
elseif(total < 346723200) % Jan 1, 1991
    delta = 6;
elseif(total < 393984000) % July 1, 1992
    delta = 7;
elseif(total < 425520000) % July 1, 1993
    delta = 8;
elseif(total < 457056000) % July 1, 1994
    delta = 9;
elseif(total < 504489600) % Jan 1, 1996
    delta = 10;
elseif(total < 551750400) % July 1, 1997
    delta = 11;
elseif(total < 599184000) % Jan 1, 1999
    delta = 12;
else
    delta = 13;
end;