function run_auto
%
% Function to run demo in automatic mode
%
global esc_flag

pause on
time_2_pause = 4;
clf

esc_flag = 0;

% Set up while loop until keystroke

% Load all .mat files up front to reduce disk reading time
% Load GLONASS az/el data
load azel_glo.mat
az_glo = az;
el_glo = el;
prn_glo = prn_gps;
mask_glo = mask;
title_azel_glo = title_str_sky;

% Load GPS az/el data
load azel_gps.mat
az_gps = az;
el_gps = el;
mask_gps = mask;
title_azel_gps = title_str_sky;

% Load GPS az/el data for a satellite
load gps_vis1.mat
sat_az_gps = az;
sat_el_gps = el;
sat_mask_gps = mask;
sat_t_los = t_los;
sat_prn_los = prn_nav_los;
sat_str_az = sprintf('GPS Azimuths for 700 km Circular Orbit at 45 Deg Inc.');
sat_str_el = sprintf('GPS Elevations for 700 km Circular Orbit at 45 Deg Inc.');
sat_str_sky = sprintf('0.5 Hrs GPS Coverage from 9/17/64 for 700 km Circular at 45 Deg Inc.');

DOP_label = {'GDOP';'PDOP';'HDOP';'VDOP';'TDOP';};

while esc_flag==0,

  %%% Satellite Visibility %%%

  % generate a GPS polar sky plot, azimuth and elevation
  % load azel_gps.mat
  if esc_flag==0,
    plot_handle = plotazel(az_gps, el_gps, prn_gps, mask_gps, title_azel_gps);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = plotaz(az_gps,el_gps,t_los_gps,prn_gps,mask_gps,title_azel_gps);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = plotel(az_gps,el_gps,t_los_gps,prn_gps,mask_gps,title_azel_gps);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    % generate a GLONASS polar sky plot, azimuth and elevation
    % load azel_glo.mat
    plot_handle = plotazel(az_glo, el_glo, prn_glo, mask_glo, title_azel_glo);
    auto_clr(time_2_pause)
  end;

  %%% Dilution of Precision for a ground stations %%%
  
  % Load the DOP plot data
  load dop_e
  
  if esc_flag==0,
    % generate the number of visible satellites
    plot(t_plot_gps,num_gps_sats,'c', ...
         t_plot_glo,num_glo_sats,'r', ...
         t_plot_all,num_all_sats,'g')
    axis([0 t_plot_gps(length(t_plot_gps)) 0 max(num_all_sats)+1])
    legend('# GPS Visible','# GLONASS Visible','# GPS & GLONASS',0);
    title(title_str_vis)
    ylabel('Number of Satellites')
    xlabel(x_str_vis)
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    % Generate the DOPS plot 
    num_plotted = 0;
    for k=1:5,
      num_plotted = num_plotted+1;
      subplot(5,1,num_plotted),
      plot(t_plot_gps,gps_dops(:,k),'c',t_plot_glo,glo_dops(:,k),'r',t_plot_all,all_dops(:,k),'g')
      axis([0 t_plot_gps(length(t_plot_gps)) 0 4])
      ylabel(DOP_label(k))
      if num_plotted == 1,
        title(title_dop_str)
      end;
    end;
    xlabel(x_dop_str)
    auto_clr(time_2_pause)
  end;

  %%% Visibility from a satellite %%%

  if esc_flag==0,
    % generate a polar sky plot from a satellite
    plot_handle = plotazel(sat_az_gps,sat_el_gps,sat_prn_los,...
                    sat_mask_gps,sat_str_sky);
    auto_clr(time_2_pause);
  end;

  if esc_flag==0,
    % generate a azimuth plot for a satellite
    plot_handle = plotaz(sat_az_gps,sat_el_gps,sat_t_los,sat_prn_los, ...
                      sat_mask_gps,sat_str_az);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    % generate a elevation plot for a satellite
    plot_handle = plotel(sat_az_gps,sat_el_gps,sat_t_los,sat_prn_los, ...
                      sat_mask_gps,sat_str_el);
    auto_clr(time_2_pause)
  end;

  %%% Dilution of Precision for a satellite %%%
  
  % Load the DOP plot data
  load dop_o
  
  if esc_flag==0,
    % generate the number of visible satellites
     plot(t_plot_gps,num_gps_sats,'c', ...
         t_plot_glo,num_glo_sats,'r', ...
         t_plot_all,num_all_sats,'g')
    % generate labels and text for the plot
    title('Number of Visible Satellites for 700 km Circular Orbit at 45 Deg Inc.')
    ylabel('Number of Satellites')
    x_text_label = sprintf('Hours Past 9/17/97 00:00:00');
    xlabel(x_text_label)
    legend('Visible GPS','Visible GLONASS','GPS & GLONASS');
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    % Generate the DOPS plot 
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
    auto_clr(time_2_pause)
  end;

  %%% Position Error Analysis/Simulation %%%

  % Show DOP contribution to position error
  if esc_flag==0,
    plot_handle = pltdoper(1);
    auto_clr(time_2_pause)
  end;

  % Show Simulated pseudo-range contribution to position error
  if esc_flag==0,
    plot_handle = pltsimpr(1);
    auto_clr(time_2_pause)
  end;

  % Show selective availability
  if esc_flag==0,
    plot_handle = pltsa(1);
    auto_clr(time_2_pause)
  end;

  % Show troposphere contribution
  if esc_flag==0,
    plot_handle = plttrop(1);
    auto_clr(time_2_pause)
  end;

  % Show ionosphere contribution
  if esc_flag==0,
    plot_handle = pltiono(1);
    auto_clr(time_2_pause)
  end;

  % Show atmospheric contribution
  if esc_flag==0,
    plot_handle = pltnav(1);
    auto_clr(time_2_pause)
  end;

  %%% Date Processing %%%
  if esc_flag==0,
    plot_handle = rdyumtxt(1);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = rdrxdata(1);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = pltnmvis(1);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = pltnmdop(1);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = pltnmpos(1);
    auto_clr(time_2_pause)
  end;

  if esc_flag==0,
    plot_handle = pltnmsa(1);
    auto_clr(time_2_pause)
  end;

  %%% Differential Processing %%%
  % Show common visibility
  if esc_flag==0,
    plot_handle = pltcovis(1);
    auto_clr(time_2_pause)
  end;

  % Pseudo range processing
  if esc_flag==0,
    plot_handle = pltnavd(1);
    auto_clr(time_2_pause)
  end;

  % Show troposphere error
  if esc_flag==0,
    plot_handle = pltrestp(1);
    auto_clr(time_2_pause)
  end;

  % Show ionosphere error
  if esc_flag==0,
    plot_handle = pltresio(1);
    auto_clr(time_2_pause)
  end;

  % DGPS Position Error Analysis
  if esc_flag==0,
    plot_handle = pltposan(1);
    auto_clr(time_2_pause)
  end;

end; % While loop





