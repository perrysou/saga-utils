function fig_handles = grnd2orb(sta_name, location, almanac_names, mask_info, ...
                                start_time, duration, time_step, plot_selection);

% fig_handles = grnd2orb(sta_name, location, almanac_names, mask_info, start_time, ...
%                           duration, time_step, plot_selection);
%
% Function to compute the visibility routines for locations fixed on 
% the Earth.  Computes the azimuth and elevation to the satellites and 
% displays an elevation vs. time, azimuth vs. time, # of satellites visible, 
% DOPS, and a sky plot showing azimuth/elevation pairs.
%
% Inputs:
%   sta_name       - Observer station names (nx?) where n is number of stations
%                    (optional), default = 'Boulder'
%   location       - Observer locations fixed on the Earth [lat, lon, hgt] (nx3)
%                    (rad,m) (optional). 
%                    Default = [40.0*pi/180, -105.16*pi/180, 1575.0] 
%   almanac_names  - Target satellites almanac names (kx?). Include the .alm extension.
%                    The file naming convention for GPS or GLONASS almanacs is 
%                    gps####.alm or glo####.alm where #### is the week number. 
%                    To specify the most recent almanac, provide only the first 3
%                    letters, just the navigation system abbreviation such as 
%                    'gps' or 'glo'. Example: ['gps'; 'glo'; 'user001';]. For user
%                    created almanacs, the full file name must be provided.
%                    (optional), default is most recent GPS almanac.
%   mask_info      - elevation masking (radians). 
%                    (1x1 to use same masking for all sites, 360 deg of azimuth).
%                    (mx4 format [sta_num elev min_az max_az] with m el/az 
%                    triples, sta_num corresponds to the sta_name input index).
%                    (mx5 format [sta_num min_elev max_elev min_az max_az] with 
%                    m el/az quadruples, sta_num corresponds to the sta_name 
%                    input index).
%                    see help on VIS_DATA for more details.
%                    (optional), default = 0 deg
%   start_time     - UTC start time [year month day hour minute second] (1x6)
%                    (optional), default = clock time of computer
%   duration       - total time duration (1x1) (sec) 
%                    (optional) default = 3600
%   time_step      - step size for computation of visibility (1x1) (sec) 
%                    (optional), default = 60
%   plot_selection - (nx9 or 1x9) optional array to limit the plots that are created
%                    each row [1=plot sky plot, 1=plot elevation, 1=plot azimuth, 
%                              1=plot # of visible satellites, 1=plot PDOP,
%                              1=plot GDOP, 1=plot HDOP, 1=plot VDOP, 1=plot TDOP]
%                    if 1x9, same options will be used for each site, 
%                    (optional), default = all plots
% Output:
%   fig_handles  - handles to the figures created
%
% See also ORB2ORB, EX_DOP_M

% Written by: Maria Evans 5/5/98
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, FIND_ALM, READYUMA, ALM2GEPH, UTC2GPS,
%                   PROPGEPH, LLA2ECEF, LOS, ECEF2NED, NED2AZEL, 
%                   VIS_DATA, NUM_VIS, PASSDATA, NED2DOPS, MAKEPLOT

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
fig_handles=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(0,8,nargin);
if ~isempty(msg)
  fprintf('%s  See help on GRND2ORB for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

if nargin < 1,  sta_name = 'Boulder'; end;
if nargin < 2,  location = [40.0*pi/180, -105.16*pi/180, 1575.0]; end;
if nargin < 3,  almanac_names = find_alm; end;
if nargin < 4,  mask_info = 0.0; end;
if nargin < 5,  start_time = fix(clock); end;
if nargin < 6,  duration = 3600; end;
if nargin < 7,  time_step = 60; end;
if nargin < 8,  plot_selection = ones(1,9); end;

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'GRND2ORB';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'sta_name';
  estruct.variable(1).req_dim = [901 size(sta_name,2);];
  estruct.variable(1).var = sta_name;
  
  estruct.variable(2).name = 'location';
  estruct.variable(2).req_dim = [901 3;];
  estruct.variable(2).var = location;
  estruct.variable(2).type = 'LLA_RAD';

  estruct.variable(3).name = 'mask_info';
  estruct.variable(3).req_dim = [1 1; 903 3; 903 4];
  if size(mask_info,2)==1,
    estruct.variable(3).type = 'ANGLE_RAD';
    estruct.variable(3).var = mask_info;
  else,
    estruct.variable(3).var = mask_info(:,2:size(mask_info,2));
    estruct.variable(3).type = 'ANGLE_RAD';
  end;

  estruct.variable(4).name = 'start_time';
  estruct.variable(4).req_dim = [1 6];
  estruct.variable(4).var = start_time;

  estruct.variable(5).name = 'duration';
  estruct.variable(5).req_dim = [1 1];
  estruct.variable(5).var = duration;
 
  estruct.variable(6).name = 'time_step';
  estruct.variable(6).req_dim = [1 1];
  estruct.variable(6).var = time_step;

  estruct.variable(7).name = 'plot_selection';
%  estruct.variable(7).req_dim = [1 9; 901 9];
  estruct.variable(7).var = plot_selection;

  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

alm_ephem = [];

% for each almanac
for i=1:size(almanac_names,1),
  alm = deblank(almanac_names(i,:));
  if size(alm,2)==3,
    [gps_alm, glo_alm] = find_alm;
    if strcmp(upper(alm),'GPS')==1,  % GPS almanac
      yuma_alm = gps_alm;
    elseif strcmp(upper(alm),'GLO')==1, % GLONASS almanac
      yuma_alm = glo_alm;
    else,
      fprintf('Invalid almanac input #%d to %s.  Returning with empty outputs.\n\n', ...
             i,estruct.func_name);
      return
    end; % if strcmp(upper(alm),'GPS')==1
  else,
    yuma_alm = alm;
  end; % if size(alm,2)==3

  alm_2_use = readyuma(yuma_alm);

  % sort out the unhealthy satellites
  I_gps_good = find(alm_2_use(:,2) == 0);
  alm_2_use = alm_2_use(I_gps_good,:);

  % convert the almanacs to ephemeris format
  [ephem] = alm2geph(alm_2_use);

  % Renumber the GLONASS satellites to 901-9??
  if strcmp(upper(yuma_alm(1:3)),'GLO')==1,  % Then GLONASS almanac
    new_num=linspace(1,100)';
    ephem(:,1)=new_num(1:size(ephem,1),1)+900*ones(size(ephem,1),1);
  end;

  % Build up the alm_ephem variables of all ephemerides
  if i==1,
    alm_ephem = ephem;
  else,
    alm_ephem(size(alm_ephem,1)+1:size(alm_ephem,1)+size(ephem,1),:)=ephem;
  end;
end;

% set start, stop, and delta-t times
start_utc = start_time;                    % set start time (UTC format)

% compute the stop time based on duration
% first convert the start time to GPS time
start_gps = utc2gps(start_utc);

% add on the duration
stop_gps(1) = start_gps(1);
stop_gps(2) = start_gps(2) + duration;

% check for week roll-overs
if stop_gps(2) > 604800
  stop_gps(1) = stop_gps(1) + fix(stop_gps(2) / 604800);
  stop_gps(2) = rem(stop_gps(2), 604800);
end % if stop_gps(2) > 604800

% compute satellite positions in ECEF 
[t_gps,prn_gps,x_gps,v_gps] = propgeph(alm_ephem,start_gps,stop_gps,time_step);

num_sta = size(location,1);

for station=1:num_sta,
 % check whether any plots are to be made
 if size(plot_selection,1) == 1,
   plot_graphs = plot_selection(1,:);
 else,
   plot_graphs = plot_selection(station,:);
 end;
 i_plot = find(plot_graphs==1);

 if any(i_plot),

  % convert the input location to radians for internal use
  here_lla(1:2) = location(station,1:2);
  here_lla(3) = location(station,3);

  % convert to ECEF
  here_ecef = lla2ecef(here_lla);

  % Get the mask info for this site and convert to radians for internal use
  clear mask
  if size(mask_info,2)==1, % Same elevation mask for all sites
    mask = mask_info;
  else,
    if size(mask_info,1)==1, % Same masking info for all sites
      mask = mask_info(1,2:size(mask_info,2));
    else, % Find correct masking info for this site
      mask_index = find(mask_info(:,1)==station);
      mask = mask_info(mask_index,2:size(mask_info,2));
    end;
  end;

  % compute LOS vectors in XYZ (ECEF)
  [t_los_gps, gps_los, los_indices, obscure_info] = ...
                       los(t_gps(1,:), here_ecef, t_gps, [prn_gps x_gps]);

  prn_gps_los = prn_gps(los_indices(:,2));

  % convert LOS in ECEF to NED
  [gps_los_ned] = ecef2ned(gps_los, here_lla);
 
  % compute az and els 
  [az el] = ned2azel(gps_los_ned);

  % call vis_data to compute visibility data
  [visible_data, I_pass] = vis_data(mask, [az, el, t_los_gps, prn_gps_los]);

  if size(visible_data,1) > 0,  % Then there is some visible data
    % compute number of visible satellites
    [t_out_vis, num_sats_vis] = num_vis(visible_data(:,3:4));

    % call passdata to generate pass numbers for each observation
    [pass_numbers] = passdata(visible_data(:,3:4), time_step*2, ...
                    [visible_data(:,5), ones(size(visible_data,1),1)*station]);

    % compute DOPS
    [gps_dops, t_gps_dops, num_gps_sats] = ned2dops(gps_los_ned(I_pass,:),visible_data(:,3:4));

    % generate a linear time scale for time plotting
    if any(t_gps_dops),
      t_dop_gps = ((t_gps_dops(:,1) -  t_gps_dops(1,1)) * 86400 * 7 + ...
                   t_gps_dops(:,2) - t_gps_dops(1,2))  / 60;
    else,
      t_dop_gps=[];
    end;
    t_plot_gps = ((t_los_gps(:,1) -  t_los_gps(1,1)) * 86400 * 7 + ...
               t_los_gps(:,2) - t_los_gps(1,2))  / 60;

    %%%% Date is computed. Now create plots.
    handles = makeplot(visible_data, pass_numbers, [t_out_vis, num_sats_vis], ...
                         [t_gps_dops, gps_dops], deblank(sta_name(station,:)), ...
                         plot_graphs, min(mask(:,1)), station-1);
    if isempty(fig_handles),
      fig_handles(1:length(handles)) = handles;
    else,
      fig_handles(length(fig_handles)+1:length(fig_handles)+length(handles)) = handles;
    end;
  end;   % if size(visible_data,1) > 0,  
 end;  % if any(iplot)
end; % End of loop over each ground station

%%%%% END ALGORITHM CODE %%%%%

% end of GRND2ORB
