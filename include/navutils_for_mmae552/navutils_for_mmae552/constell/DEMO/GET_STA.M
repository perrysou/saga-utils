function [sta_name,location,mask_info] = get_sta(sta_col)
%
% Function to convert the column strings from the
% demo into location data for the example functions
%

% Initialize all output variables
sta_name=blanks(0);
location=[];
mask_info=[];

sta_str = char(get(sta_col(1),'string'));
% Check to make sure that the station name does not exceed 12 characters
if size(sta_str,2) > 12,
   ButtonName = questdlg('Station Names are limited to 12 characters for this demo.',...
        'NOTICE','OK','Cancel');
   tmp = sta_str(:,1:12);
   tmp(size(sta_str,1)+1,1) = ' ';
   set(sta_col(1),'string',cellstr(tmp));
   sta_str = tmp(1:size(tmp,1)-1,:);
end;

lat_str = char(get(sta_col(2),'string'));
lon_str = char(get(sta_col(3),'string'));
alt_str = char(get(sta_col(4),'string'));
ele_str = char(get(sta_col(5),'string'));
azmin_str = char(get(sta_col(6),'string'));
azmax_str = char(get(sta_col(7),'string'));

rows_used=0;
for i=1:size(azmin_str,1),
  used_index = find(abs(azmin_str(i,:)) ~= 0 & abs(azmin_str(i,:)) ~= 32);
  if any(used_index), % Then an azimuth was provided for this line
    rows_used = i;
  end;
end;
   
% Initialize counters
num_sta = 0;
for i=1:rows_used,
  % Check to make sure that the sta_str even has a row for this azimuth row
  if size(sta_str,1) >= i,
    name_index = find(abs(sta_str(i,:)) ~= 0 & abs(sta_str(i,:)) ~= 32);
    if any(name_index),  % This row is a new station
      num_sta = num_sta + 1;
      sta_name(num_sta,1:length(name_index)) = sta_str(i,name_index);
      location(num_sta,:) = [str2num(lat_str(i,:)) ...
                             str2num(lon_str(i,:)) ...
                             str2num(alt_str(i,:))];             
      station_row(num_sta) = i;
    end;  % End - if any(name_index)
  end;  % End - if size(sta_str,1) >= i
end;  % End - for i=1:rows_used,
line = 0;
for i=1:num_sta,
  if i < num_sta,
    max_line = station_row(i+1);
  else,
    max_line = rows_used+1;
  end;
  for j=1:max_line-station_row(i),
    line = line + 1;
    mask_info(line,1:4) = ...
             [i str2num(ele_str(line,:)) str2num(azmin_str(line,:)),...
             str2num(azmax_str(line,:))];
  end;
end; % End - for i=1:num_sta



