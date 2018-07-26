function closewin(main_win)
%
% function to close all windows whenever a main menu option is
% changed
%
global esc_flag

open_wins = get(0,'children');
win_2_close = find(open_wins ~= main_win);
if any(win_2_close),
  close(open_wins(win_2_close));
end;

