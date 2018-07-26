function fig_handles = orb2orb(observer_almanacs, target_almanacs, mask_info, start_time, ...
                                duration, time_step, plot_selection, ...
                                obscure_altitude, earth_model_flag);

% fig_handles = orb2orb;
%        or with optional inputs
% fig_handles = orb2orb(observer_almanacs, target_almanacs, mask_info, start_time, ...
%                               duration, time_step, plot_selection, ...
%                               obscure_altitude, earth_model_flag);
%
% Function to demonstrate the visibility and DOPS routines for orbitting 
% satellites. Computes the azimuth and elevation to the GPS or GLONASS 
% satellites and displays an elevation vs. time, azimuth vs. time, # of 
% satellites visible, DOPS, and a sky plot showing azimuth/elevation pairs.
%
% Inputs:
%   observer_almanacs - Observer satellites almanac names (nx?). Include the .alm
%                       extension. Default is 'orb2orb.alm'. Optional. Can be 
%                       string array or cell array. Can also include observer 
%                       satellite names to be used only for plot titles. For example: 
%                         {'obs001.alm', 'observer_sat_1';
%                          'obs002.alm', 'observer_sat_2'}
%   target_almanacs  - Target satellites almanac names (kx?). Include the .alm extension.
%                      Default is most recent GPS almanac. Optional. Can be string arrays
%                      or cell arrays.
%                    
%   Note on almanac names: The file naming convention for GPS or GLONASS almanacs is 
%                    gps####.alm or glo####.alm where #### is the week number. 
%                    To specify the most recent almanac, provide only the first 3
%                    letters, just the navigation system abbreviation such as 
%                    'gps' or 'glo'. Example: ['gps'; 'glo'; 'user001.alm';]. For user
%                    created almanacs, the full file name must be provided.
% 
%   mask_info      - elevation masking (radians). 
%                    (1x1 to use same masking for all sites, 360 deg of azimuth).
%                    (mx4 format [sta_num elev min_az max_az] with m el/az 
%                    triples, sta_num corresponds to the sv_num input index).
%                    (mx5 format [sta_num min_elev max_elev min_az max_az] with 
%                    m el/az quadruples, sta_num corresponds to the sv_num 
%                    input index).
%                    see help on VIS_DATA for more details.
%                    (optional), default = Eliminate any observations
%                    obscured by the earth using a spherical earth model and
%                    obscure_altitude = 0.
%   start_time     - UTC start time [year month day hour minute second] (1x6)
%                    (optional), default = clock time of computer
%   duration       - total time duration (1x1) (sec) 
%                    (optional) default = 1800
%   time_step      - step size for computation of visibility (1x1) (sec) 
%                    (optional), default = 60
%   plot_selection - (nx9 or 1x9) optional array to limit the plots that are created
%                    each row [1=plot sky plot, 1=plot elevation, 1=plot azimuth, 
%                              1=plot # of visible satellites, 1=plot PDOP,
%                    1=plot GDOP, 1=plot HDOP, 1=plot VDOP, 1=plot TDOP]
%                    if 1x9, same options will be used for each observer satellite, 
%                    (optional), default = all plots
%   obscure_altitude - Altitude above earth to include in earth obscuration 
%                      (meters) (nx1). Default = 0 meters. (Optional). If not included,
%                      only the mask info will used to determine visibility. No earth 
%                      modelling will be done.
%   earth_model_flag - 0 for spherical earth, 1 for oblate earth (1x1) (default = 0)
% Output:
%   fig_handles  - handles to the figures created
%
% See also GRND2ORB, EX_DOP_M

% Written by: Maria Evans
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, READYUMA, ALM2GEPH, FIND_ALM, UTC2GPS,
%                   PROPGEPH, LOS, ECEF2NED, NED2AZEL, VIS_DATA, NUM_VIS,
%                   PASSDATA, NED2DOPS, MAKEPLOT

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
fig_handles=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(0,9,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ORB2ORB for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

if nargin < 1,  
  observer_almanacs = {'orb2orb.alm'}; 
else,
  if ~iscellstr(observer_almanacs),
     observer_almanacs = cellstr(observer_almanacs);
  end;
end;

if nargin < 2,  % No target almanacs provided
  target_almanacs = cellstr(find_alm); 
else,  % Convert string arrays to cell array format
  if ~iscellstr(target_almanacs),
     target_almanacs = cellstr(target_almanacs);
  end;
end;
if nargin < 3,  mask_info = -pi/2; end;
if nargin < 4,  start_time = fix(clock); end;
if nargin < 5,  duration = 1800; end;
if nargin < 6,  time_step = 60; end;
if nargin < 7,  plot_selection = ones(1,9); end;
if nargin < 8,  obscure_altitude = 0.0; end;
if nargin < 9,  earth_model_flag = 0; end;

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ORB2ORB';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.

  estruct.variable(1).name = 'mask_info';
  estruct.variable(1).req_dim = [1 1; 903 4; 903 5];
  estruct.variable(1).var = mask_info;

  estruct.variable(2).name = 'start_time';
  estruct.variable(2).req_dim = [1 6];
  estruct.variable(2).var = start_time;

  estruct.variable(3).name = 'duration';
  estruct.variable(3).req_dim = [1 1];
  estruct.variable(3).var = duration;
 
  estruct.variable(4).name = 'time_step';
  estruct.variable(4).req_dim = [1 1];
  estruct.variable(4).var = time_step;

  estruct.variable(5).name = 'plot_selection';
  estruct.variable(5).req_dim = [1 9; 901 9];
  estruct.variable(5).var = plot_selection;

  estruct.variable(6).name = 'obscure_altitude';
  estruct.variable(6).req_dim = [1 1];
  estruct.variable(6).var = obscure_altitude;

  estruct.variable(7).name = 'earth_model_flag';
  estruct.variable(7).req_dim = [1 1];
  estruct.variable(7).var = earth_model_flag;

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

for i=1:size(observer_almanacs,1), % For each observer almanac
  alm = char(deblank(observer_almanacs(i,1)));
  alm_2_use = readyuma(alm);
  % convert the almanacs to ephemeris format
  [ephem] = alm2geph(alm_2_use);

  % add the ephemeris to the user_ephem array
  if i==1,
    user_ephem = ephem;
  else,
    user_ephem(size(user_ephem,1)+1:size(user_ephem,1)+size(ephem,1)) = ephem;
  end;
end;

% for each target almanac
for i=1:size(target_almanacs,1),
  alm = char(deblank(target_almanacs(i,1)));
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
  else, % Renumber user constellations to 100-1??
    if strcmp(upper(yuma_alm(1:3)),'GPS')~=1,  % Not GPS or GLONASS
      new_num=linspace(1,size(ephem,1),size(ephem,1))';
      ephem(:,1)=new_num(1:size(ephem,1),1)+100*ones(size(ephem,1),1);
    end;
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

% compute user satellite positions in ECEF 
[t_user,prn_user,x_user,v_user] = ...
         propgeph(user_ephem,start_gps,stop_gps,time_step);

% compute navigation satellite positions in ECEF 
[t_gps,prn_gps,x_gps,v_gps] = propgeph(alm_ephem,start_gps,stop_gps,time_step);

DOP_label = ['GDOP';'PDOP';'HDOP';'VDOP';'TDOP';];
n_plots = 0;
   
% compute LOS vectors in XYZ (ECEF) for this satellite to the GPS constellation
[t_los_all, los, los_indices, obscure_info] = ...
       los(t_user, [prn_user x_user], t_gps, [prn_gps x_gps], earth_model_flag);

% Get arrays from los output indices
prn_user_los = prn_user(los_indices(:,1));
user_pos_out = x_user(los_indices(:,1),:);
prn_nav_los = prn_gps(los_indices(:,2));
nav_pos_out = x_gps(los_indices(:,2),:);

% convert LOS in ECEF to NED   
[los_ned] = ecef2ned(los, ecef2lla(user_pos_out));

% compute az and els 
[azimuth elevation] = ned2azel(los_ned);

% Loop over the user input satellites for plotting
for k=1:size(user_ephem,1),  

  % get the indices that match this satellite
  plot_prn = find(prn_user_los == user_ephem(k,1));
  az = azimuth(plot_prn);
  el = elevation (plot_prn);
  t_los = t_los_all(plot_prn,:);
  prn_nav = prn_nav_los(plot_prn);
  prn_user = prn_user_los(plot_prn);
  ned_los = los_ned(plot_prn,:);

  % Get the mask info for this satellite and convert to radians for internal use
  clear mask
  if size(mask_info,2)==1, % Same elevation mask for all satellites
    mask = mask_info;
  else,
    if size(mask_info,1)==1, % Same masking info for all satellites
      mask = mask_info(1,2:size(mask_info,2));
    else, % Find correct masking info for this site
      mask_index = find(mask_info(:,1)==station);
      mask = mask_info(mask_index,2:size(mask_info,2));
    end;
  end;

  % call vis_data to compute visibility data
  if nargin < 8 & nargin >= 3,  % Use only mask data and don't obscure out the Earth
    [visible_data, I_pass] = vis_data(mask, [az, el, t_los, prn_user, prn_nav]);
  else,  % Then mask out Earth from visibility
    [visible_data, I_pass] = vis_data(mask, [az, el, t_los, prn_user, prn_nav],...
                   obscure_info(plot_prn,:), obscure_altitude);
  end;
  
  if size(visible_data,1) > 0,  % Then there is some visible data
    % compute number of visible satellites
    [t_out_vis, num_sats_vis, id_num_vis] = num_vis(visible_data(:,3:4), visible_data(:,5));

    % call passdata to generate pass numbers for each observation
    [pass_numbers, pass_times] = passdata(visible_data(:,3:4), time_step*2, ...
                     [visible_data(:,5), visible_data(:,6)]);

    % compute DOPS
    [gps_dops, t_gps_dops, num_gps_sats] = ned2dops(ned_los(I_pass,:), ...
                                          t_los(I_pass,:));

    if any(plot_prn), % there is something to plot
      % check whether any plots are to be made
      if size(plot_selection,1) == 1,
        plot_graphs = plot_selection(1,:);
      else,
        plot_graphs = plot_selection(station,:);
      end;

      % generate a linear time scale for time plotting
      if any(t_gps_dops),
        t_dop_gps = ((t_gps_dops(:,1) -  t_gps_dops(1,1)) * 86400 * 7 + ...
                   t_gps_dops(:,2) - t_gps_dops(1,2))  / 60;
      else,
        t_dop_gps=[];
      end;
      t_plot = ((t_los(:,1) -  t_los(1,1)) * 86400 * 7 + ...
                   t_los(:,2) - t_los(1,2))  / 60;

      %%%% Data is computed. Now create plots.
      if size(observer_almanacs,2)==2,
        label_string = sprintf('%s',char(observer_almanacs(k,2)));
      else,
        label_string = sprintf('User Satellite %d',user_ephem(k,1));
      end;
      handles = makeplot([visible_data(:,1:4), visible_data(:,6)], pass_numbers, ...
                           [t_out_vis, num_sats_vis], ...
                           [t_gps_dops, gps_dops], label_string, ...
                           plot_graphs, min(visible_data(:,2)), k-1);
      if isempty(fig_handles),
        fig_handles(1:length(handles)) = handles;
      else,
        fig_handles(length(fig_handles)+1:length(fig_handles)+length(handles)) = handles;
      end;
    end;  % end if any(plot_prn)
  end;  %   if size(visible_data,1) > 0, 
end;  % end for k=1:size(user_ephem,1)

% end of ORB2ORB
