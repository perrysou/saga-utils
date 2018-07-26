function demo_orb(flag,main_win)
global vis_handles show_flag
%
% Call ex_vis_o or ex_dop_o to demonstrate the visibility or DOP
% from an orbiting satellite
%
% Input:
%    flag = 1 for Visibility plots
%              = 2 for DOP plots
%    main_win = optional input included only for manual demo mode
%

show_flag = flag;

clear vis_handles

%%% The following code was run initially to create gps_vis1.mat
%%% and glo_vis1.mat. These are loaded and plotted for the demo
%user_orbit = [1 7000000 0 45 0 0 0 923 0];
%user_mask = -25;
%nav_system = 1;  % GPS
%start_time = [1997 9 17 0 0 0];
%vis_handles = ex_vis_o(user_orbit, user_mask, nav_system, start_time);
%%%

% DEMO Version

% Get the screen size in pixels to use for location of plots
set(0,'units','pixels');
screen_size = get(0,'screensize');
y_max = screen_size(2) + screen_size(4) - 60;
x_max = screen_size(1) + screen_size(3) - 50;
x_step = 110;
y_step = 60 + y_max/2;

if flag==1,
  % Load GPS visibility data
  load gps_vis1
else,
  % Load the DOPS data
  load dop_o
end;

% Set the main window to visible again
if nargin == 2,
  set(main_win,'visible','on');
end;

n_plots = 0;

if flag==1,

  % generate the plot of azimuth versus time                    
  [azel_vis, I_pass] = vis_data(0,[az el]);
  id_nums = [ones(size(prn_nav_los(I_pass))) prn_nav_los(I_pass)];
  [pass_numbers, pass_times] = passdata(t_los(I_pass,:), 600,id_nums); 
  [t_vis, num_sats] = num_vis(t_los(I_pass,:));
  visible_data = [azel_vis t_los(I_pass,:) prn_nav_los(I_pass)];   
  dops = zeros(10,7);
  plot_selection = [1 1 1 1 1 1 1 1 1];
  plot_handle = makeplot(visible_data,pass_numbers,[t_vis num_sats],...
                         dops, 'Example Orbit to GPS Satellite VIsibility',...
                         plot_selection);

else,  % Show DOP plots
  x_step = x_max/2;
  y_step = 60 + y_max;
  x_min = x_max/2;
  y_min = 30;
  DOP_label = {'GDOP';'PDOP';'HDOP';'VDOP';'TDOP';};

  n_plots = n_plots + 1;
  vis_handles(n_plots) = figure('position',[x_min-x_step y_min x_max/2 y_max/2], ...
       'NumberTitle','off','Units','Pixels','Color','Black');
  colordef none
  set(vis_handles(n_plots),'Name',sprintf('# Visible Satellites for 700 km Circular Orbit at 45 Deg Inc.'));
  plot(t_plot_gps,num_gps_sats,'c',t_plot_glo,num_glo_sats,'r',t_plot_all,num_all_sats,'g')

  % generate labels and text for the plot
  title('Number of Visible Satellites for 700 km Circular Orbit at 45 Deg Inc.')
  ylabel('Number of Satellites')
  x_text_label = sprintf('Hours Past 9/17/97 00:00:00');
  xlabel(x_text_label)
  legend('Visible GPS','Visible GLONASS','GPS & GLONASS');

  % generate the DOP plots
  n_plots = n_plots + 1;
  fig_handles(n_plots) = figure('position',[x_min y_min x_max/2 y_max/1.2], ...
     'NumberTitle','off','Units','Pixels','Color',[0 0 0],...
     'Name','DOPs');
  colordef none

  num_plotted = 0;
  for k=1:5,
    num_plotted = num_plotted+1;
    subplot(5,1,num_plotted),
    plot(t_plot_gps,gps_dops(:,k),'c',t_plot_glo,glo_dops(:,k),'r',t_plot_all,all_dops(:,k),'g')
    axis([0 t_plot_gps(length(t_plot_gps)) 0 6])
    ylabel(DOP_label(k))
    if num_plotted == 1,
      title_dop_str = sprintf('700 km Circular, 45 Deg Inc., Blue=GPS, Red=GLONASS, Green=Both');
      title(title_dop_str)
    end;
  end;
  x_dop_str = sprintf('Hours Past 9/17/97 00:00:00');
  xlabel(x_dop_str)

end;

