function del_sta(sta_col,sta_2_del)
%
% Function to add a station to the table or add more masking
% to an existing station
% Input:
%   sta_col   = uicontrol for the columns of the input table
%   sta_2_del = Number of station to delete from list
%
[sta_name,location,mask_info] = get_sta(sta_col);
num_sta = size(sta_name,1);
num_line = size(mask_info,1);

% Determine the line numbers to delete
if length(mask_info) > 0,
 i_match = find(mask_info(:,1)==sta_2_del);
 if any(i_match),
  for j=1:7,
    old_str = get(sta_col(j),'string');
    if i_match(1) == 1,
      first_line = i_match(length(i_match))+1;
      new_str = old_str(first_line:size(mask_info,1),:);
    else,
      new_str = old_str(1:i_match(1)-1,:);
      if i_match(length(i_match)) < num_line;
        new_str(size(new_str,1)+1:num_line-length(i_match),:) = ...
            old_str(i_match(length(i_match))+1:size(mask_info,1),:);
      end;
    end;

    set(sta_col(j),'string',blanks(0));  % First blank out the entire string
    set(sta_col(j),'string',new_str);
  end;
 end;
end;
