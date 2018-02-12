% Function hhmm = dechr2hhmm(hour) converts a decimal hour of day into a
% string hhmm, expressing the time in hours and minutes.
%
% Input:
%   hour    - decimal valued hour, e.g., 5.5.
% Output:
%   hhmm    - string output of time, e.g., '0530'.
%
% Seebany Datta-Barua
% 30 Apr 2008
% 1 May 2008 Vectorize to accept array input.
% 18 Mar 2008 Take care of not quite round input hours.

function hhmm = dechr2hhmm(hour)

% Create the string of the filename we need to open.
hh = floor(hour);
% Convert fractional hour to minutes and make string.
mm = round(rem(hour, max(hh,1))*60);

% Round if mm is 60.
rows = find(mm == 60);
hh(rows) = hh(rows)+1;
mm(rows) = 0;

hhstr = num2str(hh);
rows = find(hh < 10);
hhstr(rows,1:2) = [repmat('0',size(hhstr(rows,1),1),1) hhstr(rows,1)];

mmstr = num2str(mm);
rows = find(mm < 10);
mmstr(rows,1:2) = [repmat('0',size(mmstr(rows,1),1),1) mmstr(rows,1)];

hhmm = [hhstr mmstr];