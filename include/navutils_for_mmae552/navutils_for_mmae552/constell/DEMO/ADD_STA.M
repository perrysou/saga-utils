function add_sta(sta_col,sta_2_mod)
%
% Function to add a station to the table or add more masking
% to an existing station
% Input:
%   sta_col   = uicontrol for the columns of the input table
%   sta_2_mod = optional, include only to add more masking
%               Set equal to the station to modify
%
[sta_name,location,mask_info] = get_sta(sta_col);
num_sta = size(sta_name,1);
num_line = size(mask_info,1);
if num_line==11,  % Limit of 11 lines for demo
  ButtonName = questdlg('Limited to 11 lines of data for Demo.', ...
                        'NOTICE','OK','Cancel');
  return;
end;

if nargin==1,
  default=setstr(ones(7,30)*32);
  default(1,1:8) = 'New Site';
  default(2,1:4) = '0.00';
  default(3,1:4) = '0.00';
  default(4,1:3) = '0.0';
  default(5,1:3) = '5.0';
  default(6,1:3) = '0.0';
  default(7,1:5) = '360.0';

  for j=1:7,
    new_str = char(get(sta_col(j),'String'));
    add_str = deblank(default(j,:));
    new_str(num_line+1,1:length(add_str)) = add_str;
    set(sta_col(j),'String',cellstr(new_str));
  end;
else,
  default=setstr(ones(7,30)*32);
  default(5,1:3) = '5.0';
  default(6,1:3) = '0.0';
  default(7,1:5) = '360.0';
  if length(mask_info) > 0,
   i_match = find(mask_info(:,1)==sta_2_mod);
   if any(i_match),
    line_2_mod = i_match(length(i_match))+1;

    for j=1:7,
      new_str = char(get(sta_col(j),'String'));
      insert_str = deblank(default(j,:));
      new_str(line_2_mod+1:num_line+1,:) = ...
              new_str(line_2_mod:num_line,:);
      add_str = default(j,:);
      nchar = min(length(add_str),size(new_str,2));
      if any(add_str),
        new_str(line_2_mod,1:nchar) = add_str(1:nchar);
      end;
      set(sta_col(j),'String',cellstr(new_str));
    end;
   end;
  end;
end;

