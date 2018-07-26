function auto_clr(time_2_pause)
%
% function to pause, clear the current figure, then reset to black
% Input:
%   time_2_pause = number of seconds to pause
%
global esc_flag

for i=1:time_2_pause,
  if esc_flag==0,
    pause(1)
  end;
end;
clf
colordef(gcf,'black')
set(gcf,'color','black')


