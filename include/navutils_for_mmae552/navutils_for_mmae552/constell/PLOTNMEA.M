function [fig_han,gga_out,gsv_out,gsa_out] = ...
                      plotnmea(file_name,plot_flags,edit_data_flag);

% [fig_han,gga_data,gsv_data,gsa_data] = plotnmea(file_name,plot_flags,edit_flag);
%
% Example program for reading and plotting NMEA data.               
%
% Inputs:
%   file_name  - name of file to read NMEA data from (must be a string) (optional)
%                 default = 'nmeadata.dat' (provided with the Toolbox).
%   plot_flags - 1x3 vector indicating which plots to generate (optional)
%                [plot_gga_data plot_gsv_data plot_gsa_data],
%                1 = plot the data, 0 = do not generate plots for that message
%                default = [1 1 1]
%   edit_flag  - flag indicating if the NMEA data is to be edited before
%                 plotting.  (optional) 1 = edit, 0 = do not edit, 
%                 default = 1.  Editing is based on the lat/lon not set to
%                 zero values, HDOP < 15, and 4 or more satellites visible.
% Outputs: 
%   fig_han    - handles to the figures created by the example function
%   gga_out    - NMEA data from $GPGGA messages
%                 gga_out is an nx12 matrix with columns
%                [hours, mins, secs, lat, lon, hgt, fix_qual, num_svs, HDOP, ...
%                 geoid_height, DGPS_latency, DGPS_station_number]
%                lat and lon are in deg, hgt and geoid_height are in m
%                fix_quality - 1 = GPS, 2 = DGPS
%
%   gsa_out    - NMEA data from $GPGSA messages
%                 gsa_out is an nx17 matrix with columns
%                 [auto, fix_dim, prn_1, prn_2, ..., prn_12, PDOP, HDOP, VDOP]
%                 auto - 0 = auto 2/3-D mode, 1 = manual 2/3-D mode
%                 fix_dim is the dimension fix for this data, 2 = 2-D, 3 = 3-D 
%         
%   gsv_out    - NMEA data from $GPGSV messages
%                 gsv_out is an nx49 matrix with columns
%                 [num_sats, prn_1, el_1, az_1, snr_1, prn_2, ..., snr_12]
%                 num_sats - number of satellite being tracked
%                 elevation (el) and azimuth (az) are in deg 
%                 signal-to-noise ration in dB (definition may vary by receiver)
%
% See also READNMEA, PARSNMEA, PARSEGSA, PARSEGSV 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: READNMEA, PLOTSKY, LLA2ECEF, ECEF2NED

if nargin < 1
  file_name = 'nmeadata.dat';
end % if nargin < 1

if nargin < 2
  plot_flags = [1 1 1];
end % if nargin < 2

if nargin < 3
  edit_flag = 1;
end % if nargin < 3

% set plotting options
% flag to generate GGA message plots (1 = generate)
generate_gga_plots = plot_flags(1);     

% flag to generate GSV message plots (1 = generate)
generate_gsv_plots = plot_flags(2);    

% flag to generate GSA message plots (1 = generate)
generate_gsa_plots = plot_flags(3);    

% set the data editing flag
edit_data = edit_flag;      % 1 = data editing, 0 = plot all data

% read in NMEA data
[gga_out, gsa_out, gsv_out] = readnmea(file_name);

% set the figure counter to 1
fig_count = 1;
fig_han = [];

% edit the data for bad positions if desired
if edit_data == 1 & ~isempty(gga_out)                            
  % select based on lat and lon not equal to zero (a defualt value)
  % and at least 4 satellites in use and HDOP < 15
  I_good = find(gga_out(:,4) ~= 0 & gga_out(:,5) ~= 0 & ...
                gga_out(:,8) > 3 & gga_out(:,9) < 15);
  gga_out_all = gga_out;
  gga_out = gga_out(I_good,:); 
end % if edit_data == 1

if generate_gsv_plots == 1 & ~isempty(gsv_out)
  % format the prn/az/el data for the plotsky function
  % search over the satellite number fields in the GSV messages 
  count = 1;                   % start index
  for i = 1:size(gsv_out,1)    % loop over each gsv message
    for j = 1:gsv_out(i,1)     % loop over the number of satellites in each message
      sv_index = (j - 1) * 4 + 2;
      prn(count) = gsv_out(i,sv_index + 0);   % prn
      el(count) = gsv_out(i,sv_index + 1);    % el
      az(count) = gsv_out(i,sv_index + 2);    % az
      snr(count) = gsv_out(i,sv_index + 3);   % snr  
      count = count + 1;
    end % for j = 1:gsv_out(i,1)
  end % for i = 1:size(gsv_out,1)

  % conver az and el to radians
  d2r = pi / 180;
  az = az * d2r;
  el = el * d2r;
    
  % get the matrices in the right format (nx1 not 1xn)
  el = el';
  az = az';
  prn = prn';
  snr = snr';
  
  % now create a figure
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;
  
  % get times from the gga messages because GSA does not have time 
  hours = gga_out(:,1) + gga_out(:,2) / 60 + gga_out(:,3) / 3600;  % plot time
  
  % truncate gsv_out messages to be the same size as the hours (gga_out)
  lg = length(hours);
  
  % verify that the length of the gga and gsv messages are the same
  % if data is missed or a gga message is the last in the file this could happen
  % if data has been skipped, warn the user that there could be a time tagging
  % error
  if size(gsv_out,1) ~= lg
    if size(gsv_out,1) ~= lg - 1 
      fprintf('Warning from PLOTNMEA.  There seems to be data missing from the\n');
      fprintf('NMEA data file %s.  There should be a GGA message for \n',...
               file_name);
      fprintf('every GSV message.  Check the data file for missing data.\n');
      fprintf('%d GGA messages were found and %d GSV messages were found.\n',...
               size(gga_out,1),size(gsv_out,1));
    end % if size(gsv_out,1) ~= lg - 1 
    
    % reduce the size of lg to be the same as the GSV messages
    
    lg = size(gsv_out,1);
    
    if lg > size(hours,1)
      lg = size(hours,1);
    end % if lg > size(hours,1)
    
    hours = hours(1:lg);
    hours_nvis = hours;
    
  end % if size(gsv_out,1) ~= lg
  
  nvis = gsv_out(1:lg,1); 
  plot(hours_nvis,nvis);
  
  % set the axis to show the number of satellites (add +/-1 to the y-axis)
  vs = axis;
  axis([vs(1) vs(2) vs(3)-1 vs(4)+1]);
  
  title('Number of Satellites in View (from GPGSV Message)')
  xlabel('Hours past start of the day (UTC)')
  ylabel('Number of GPS Satellites Visible')
  set(gcf,'Name','Number of GPS Satellites Visible')

  % generate a sky view plot of the data
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  % set the title string for the sky plot
  title_string = ...
   sprintf(...
   'Sky Plot from GSV Data Starting at %4.2f (hrs) with %4.2f (hrs) Duration',...
          hours(1),hours(lg) - hours(1)); 
  % set the masking for the sky plot
  sky_mask = -pi/2;
                            
  plotsky(az,el,prn,title_string);  
  
  clear hours lg
  
end % if generate_gsv_plots == 1

if generate_gsa_plots == 1  & ~isempty(gsa_out)
  % get times from the gga messages because GSA does not have time 
  hours = gga_out(:,1) + gga_out(:,2) / 60 + gga_out(:,3) / 3600;  % plot time
  
  lg = size(hours,1);
  
  % verify that the length of the gga and gsv messages are the same
  % if data is missed or a gga message is the last in the file this could happen
  % if data has been skipped, warn the user that there could be a time tagging
  % error
  if size(gsa_out,1) ~= lg
    if size(gsa_out,1) ~= lg - 1 
      fprintf('Warning from PLOTNMEA.  There seems to be data missing from the\n');
      fprintf('NMEA data file %s.  There should be a GGA message for \n',...
               file_name); 
      fprintf('every GSA message.  Check the data file for missing data.\n');
      fprintf('%d GGA messages were found and %d GSA messages were found.\n',...
               size(gga_out,1),size(gsa_out,1));
    end % if size(gsa_out,1) ~= lg - 1 
    
    % reduce the size of lg to be the same as the GSV messages
    
    lg = size(gsa_out,1);

    if lg > size(hours,1)
      lg = size(hours,1);
    end % if lg > size(hours,1)

    hours = hours(1:lg);
  end % if size(gsv_out,1) ~= lg
  
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;
  
  set(gcf,'Name','GPGSA Message DOPs')

  subplot(3,1,1)
  plot(hours,gsa_out(1:lg,15))
  ylabel('PDOP')
  title('DOPs Read from GPGSA Messages')

  subplot(3,1,2)
  plot(hours,gsa_out(1:lg,16))
  ylabel('HDOP')

  subplot(3,1,3)
  plot(hours,gsa_out(1:lg,17))
  ylabel('VDOP') 
  xlabel('Hours past start of the day (UTC)') 
  
  hours_dop = hours;
  pdop = gsa_out(1:lg,15);
  hdop = gsa_out(1:lg,16);
  vdop = gsa_out(1:lg,17);
  
end % if generate_gsa_plots == 1  

if generate_gga_plots == 1  & ~isempty(gga_out)

  hours = gga_out(:,1) + gga_out(:,2) / 60 + gga_out(:,3) / 3600;  % plot time
  hours_pos = hours;

  % generate a plot in 3-D of the positions over that time  
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  fig1 = gcf;
  plot3(gga_out(:,5),gga_out(:,4),gga_out(:,6))
  grid
  xlabel('longitude')
  ylabel('latitude')
  zlabel('height (m)')
  title('3D Display of NMEA Logged Positions')
  set(gcf,'Name','3-D Logged Positions')

  % now look at the data in NED coordinates
  % get the lat, lon, hgt in rad  
  clear lla dx
  lla = zeros(size(gga_out,1),3);
  lla(:,1:2) = gga_out(:,4:5) * pi / 180;
  lla(:,3) = gga_out(:,6);

  x_ecef = lla2ecef(lla);   % convert to ECEF

  mean_lla = mean(lla);     % compute the mean lla
  mean_ecef = lla2ecef(mean_lla);  % mean in ECEF
 
  dx(:,1) = x_ecef(:,1) - mean_ecef(1);   % x-delta
  dx(:,2) = x_ecef(:,2) - mean_ecef(2);   % y-delta
  dx(:,3) = x_ecef(:,3) - mean_ecef(3);   % z-delta

  % convert to NED
  d_ned = ecef2ned(dx,mean_lla); 

  % generate plot
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  fig2 = gcf;
  plot3(d_ned(:,1),d_ned(:,2),d_ned(:,3))
  grid               
  xlabel('delta-N (m)')
  ylabel('delta-E (m)')
  zlabel('delta-D (m)')
  title('3D Display of Position Variation')
  set(gcf,'Name','3-D Position Variation from the Mean')

  % offset figure 2
  x_pos = get(fig1,'Position');
  set(fig2,'Position',[x_pos(1)+20 x_pos(2)-20 x_pos(3:4)]) 

  % generate a figure with the number of satellite 
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  fig3 = gcf;
  hours = gga_out(:,1) + gga_out(:,2) / 60 + gga_out(:,3) / 3600;  % plot time
  plot(hours, gga_out(:,8)) 
  xlabel('Hours past the start of the day')
  ylabel('Number of Satellites')
  title('Number of GPS Satellites Tracked')
  set(gcf,'Name','Number of GPS Satellites Tracked')
  
  % set the axis size to be +/- 1 from what it is by defaultfor easier viewing
  va = axis;
  axis([va(1) va(2) va(3)-1 va(4)+1]);
  
  % offset figure 3
  x_pos = get(fig2,'Position');
  set(fig3,'Position',[x_pos(1)+20 x_pos(2)-20 x_pos(3:4)]) 

  % generate a figure with the HDOP
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  fig4 = gcf;
  plot(hours, gga_out(:,9)) 
  xlabel('Hours past the start of the day')
  ylabel('HDOP')
  title('HDOP')
  set(gcf,'Name','HDOP from GGA Messages')

  % offset figure 4
  x_pos = get(fig3,'Position');
  set(fig4,'Position',[x_pos(1)+20 x_pos(2)-20 x_pos(3:4)]) 

  % compute sigmas
  sigmas = std(d_ned); 

  % figure 5 Horizontal position
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  plot(gga_out(:,5),gga_out(:,4))
  xlabel('Longitude (degrees)')
  ylabel('Latitude (degrees)') 
  title('Horizontal Position Solutions')

  set(gcf,'Name','Horizontal Position Solutions')

  % figure 6 Horizontal position error
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  plot(d_ned(:,2),d_ned(:,1))
  xlabel('East Variation (m)')
  ylabel('North Variation (m)') 
  title('Horizontal Position Variations')
  sigma_n_text = sprintf('Sigma North = %8.2f',sigmas(1));
  sigma_e_text = sprintf('Sigma East = %8.2f',sigmas(2));
  legend('Lat/Lon Positions',sigma_n_text,sigma_e_text,0)

  set(gcf,'Name','Horizontal Position Solutions')

  % figure 7, North, East and Down deviations as a function of time
  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  set(gcf,'Name','NED Deviations')

  subplot(3,1,1)
  plot(hours,d_ned(:,1))
  ylabel('North (m)')
  title('NED Deviations from the Mean')

  subplot(3,1,2)
  plot(hours,d_ned(:,2))
  ylabel('East (m)')

  subplot(3,1,3)
  plot(hours,d_ned(:,3))
  ylabel('Down (m)')
  xlabel('Hours past 0:00 Feb 1 (UTC)')

  % figure(8), Power Spectral Densities
  num_points = size(d_ned,1);
  num_freqs = fix(num_points / 2) - 1;    % don't show the DC component (f = 0)

  Y = fft(d_ned);      % 1-D fft on each column

  Pyy = Y .* conj(Y) / num_points;    % PSD

  f = 1 * (1:num_freqs) / num_points;    % use the first half with 1/sec data
  period = 1 ./ f;    % plot period instead of frequency
  period = period / 60;    % convert to minutes (from seconds)

  fig_han(fig_count) = figure;
  fig_count = fig_count + 1;

  set(gcf,'Name','NED Error Power Spectral Densities')

  subplot(3,1,1)
  plot(period,Pyy(1:num_freqs,1)) 
  ylabel('Power (North)')
  title('NED Variation Power Spectral Densities')

  subplot(3,1,2)
  plot(period,Pyy(1:num_freqs,2)) 
  ylabel('Power (East)')

  subplot(3,1,3)
  plot(period,Pyy(1:num_freqs,3)) 
  ylabel('Power (Down)')
  xlabel('Period (minutes)')
end % if generate_gga_plots == 1

save exnmedat

% decimate the data if storing output for the demo with long time spans
use_for_demo = 1;
if use_for_demo == 1
  data_freq = 10;
  data_freq_power = 10;     % use if the data is already decimated
  hours_dop = hours_dop(1:data_freq:length(hours_dop));
  pdop = pdop(1:data_freq:length(pdop)); 
  hdop = hdop(1:data_freq:length(hdop)); 
  vdop = vdop(1:data_freq:length(vdop)); 
  hours_nvis = hours_nvis(1:data_freq:length(hours_nvis));
  nvis = nvis(1:data_freq:length(nvis)); 
  hours_pos = hours_pos(1:data_freq:length(hours_pos)); 
  d_ned = d_ned(1:data_freq:length(d_ned),:); 
  
  % recompute the PSD stuff
  clear Pyy period num_freqs
                             
  % figure(8), Power Spectral Densities
  num_points = size(d_ned,1);
  num_freqs = fix(num_points / 2) - 1;    % don't show the DC component (f = 0)

  Y = fft(d_ned);      % 1-D fft on each column

  Pyy = Y .* conj(Y) / num_points;    % PSD

   % use the first half with 1/data_freq data
  f = 1 * (1:num_freqs) / num_points / data_freq_power;   
  period = 1 ./ f;    % plot period instead of frequency
  period = period / 60;    % convert to minutes (from seconds)

  % save the demo data file
  save demonmea hours_nvis nvis hours_dop pdop hdop vdop hours_pos d_ned period Pyy num_freqs
end % if use_for_demo == 1  

% end of PLOTNMEA

