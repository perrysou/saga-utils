function [t_pr,prn,pr,pr_errors,pr_ndex,obscure_info]=...
               pseudo_r(t_user,x_user,v_user,t_gps,x_gps,v_gps,model,...
                        seed,code_noise,carrier_noise,ephem,model_data)

% [t_pr,prn,pr] = pseudo_r(t_user,x_user,v_user,t_gps,x_gps,v_gps);
%                      or
% [t_pr,prn,pr,pr_errors,obscure_info] = ...
%   pseudo_r(t_user,x_user,v_user,t_gps,x_gps,v_gpsmodel,seed,code_noise, ...
%            carrier_noise,ephem,model_data);
%
% Computes pseduo-range (pr) and accumulated carrier phase (cph) measurements
% given a user trajectory and the GPS/GLONASS satellite positions. 
% Errors that can be included in the measurements include S/A (epislon and 
% dither), troposphere, ionosphere, receiver clock bias and clock drift,  
% receiver code and carrier tracking noise, satellite motion, GPS satellite
% clock, and relativity.  Pseudorange and phase data are generated for all
% satellites and user positions.  Editing and masking of the measurement
% data can be accomplished using the function VIS_DATA.
%
% Input:
%   t_user - GPS time vector for user trajectory [GPS_week, GPS_sec] (nx2)
%             There must be coincide times with the times in the GPS time vector. 
%             If there are no coincident times, no pseudo-range or phase data
%             will be generated and the output matrices will be empty.
%             valid GPS_week values are 1-3640 (years 1980-2050)
%             valid GPS_sec values are 0-604799
%   x_user - ECEF/ECI position vectors for the user vehicle [x,y,z] (nx3) (m) or
%             if more than 1 receiver is being used, [receiver_num x y z] (nx4) or
%             if multiple antenna are used, [receiver_num ant_num x y z] (nx5)
%   v_user - ECEF/ECI velocity vectors for the user vehicle [x,y,z] (nx3) (m/s)
%   t_gps  - GPS time vector for GPS positions [GPS_week, GPS_sec] (mx2)
%             valid GPS_week values are 1-3640 (years 1980-2050)
%             valid GPS_sec values are 0-604799
%   x_gps  - ECEF/ECI position vectors for GPS satellites [prn,x,y,z] (mx4) (m)
%   v_gps  - ECEF/ECI velocity vectors for GPS satellites [prn,x,y,z] (mx3) (m/s)
%   model  - flags controlling which contributions to the 
%             PR errors are modeled (optional).  (1x11)
%              [sa_eps dither troposphere ionosphere receiver_clock 
%               receiver_noise line_bias sat_motion sat_clock earth_rotation
%               relativity]
%              a value of 1 indicates useage of the model
%              and a value of zero indicates no use of that model.  Use values
%              of 2 to implement user supplied models.  See the code for where
%              to insert the user models.  A warning is given if a user model
%              is selected and none is supplied.
%              Default = [1 1 1 1 1 1 1 0 0 0 0 0].
%   seed   - seed value for random number generator (optional)
%             Default value is 0. 
%   code_noise    - 1 sigma estimate of the receiver code noise (optional) (1x1)
%                    (m), default = 1
%   carrier_noise - 1 sigma estimate of the receiver carrier noise (optional) 
%                    (1x1) (m), default = 0.01  
%   ephem         - ephemeris matrix for all satellites (nx24). (optional)
%                    Used to compute satellite clock. If not provided, no 
%                    GPS satellite clock effects will be computed
%                    The columns of ephemeris matrix are ...
%                    [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                    i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                    GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                    Ephemeris parameters are from ICD-GPS-200 with the 
%                    exception of perigee_rate. 
%   model_data    - structure with data for each of the models.  (optional)
%                   If not provided, model defaults will be used.  Valid model
%                   structure elements are ...
%                    model_data.sa_eps (1x1) see SA_EPS for details
%                    model_data.sa_dither (1x3) see SA_CLOCK for details
%                    model_data.tropo (1x3) see TROPDLAY for details
%                    model_data.iono (1x8) see IONODLAY for details
%                    model_data.rec_clock (1x2) see CLOCKERR for details
%                    model_data.line_bias (nx3) [rec_num ant_num line_bias_sigma] 
%                                default = 1 meter for all line biases
%
% Note:  If a single user time and position set is input with multiple GPS
%        satellite time/positions, the user is considered to be stationary and
%        PR and CPH measurements will be generated for all GPS times.  For
%        this special case, no obscure_info will be returned (obscure_info = []).
%
% Output:
%   t_pr   - GPS time associated with the PR, [GPS_week GPS_sec] (kx2)
%              k = num_time_steps x number of visible satellites
%   prn    - satellite number for this pseudo-range (kx1)
%   pr     - pseudo-range and accumulated carrier phase measurements associated 
%             with the corresponding t_pr and prn (kx2) [pr cph] (m and cycles)
%   pr_errors - modeled errors added to the geometric range to obtain
%                a pseudorange and carrier phase measurement (m or m/s) (kx13)
%                [sa_eps_err sa_clk_err trop_wet trop_dry iono ...
%                 clk_bias clk_drift code_white carrier_white ...
%                 line_bias sat_motion sat_clock relative]    
%   pr_ndex   - indices to obtain relationship between output PR matrix and
%                the input positions and velocity [x_user_ind, x_gps_ind] (kx2)
%   obscure_info - contains information needed to determine whether the earth
%                  obscures the line-of-sight (kx3) [tangent_altitude, alpha, beta].
%                  Alpha is the angle between the x1 and x2 vectors. 
%                  Beta is the angle between the x1 vectors and the radius to the
%                  tangent point of the los vectors.
%                  An observation is obscured if the tangent vector magnitude is
%                  below the users tolerance, and alpha is greater than beta. 
%
% See also DOPPLER, SA_EPS, SA_CLOCK, CLOCKERR, TROPDLAY, IONODLAY, VIS_DATA

% References: 
% ICD-GPS-200C-002, 25 Sept, 1997, pages 88-89, 102.
% 
% Global Positioning System: Theory and Applications
% Volume 2, Parkinson and Spilker, pages 8, 28.
%
% 'GPS Theory and Practice", Hoffman-Wellenhoff, pages 89-90.

% Written by: Jimmy LaMance 1/15/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, LOS, NORMVECT, SA_EPS, SA_CLOCK, ADDAMBIG,
%                   CLOCKERR, ECEF2LLA, ECEF2NED, NED2AZEL, TROPDLAY, 
%                   IONODLAY, GPST2SEC, KEPLR_EQ

% WGS-84 constants
LIGHT_SPEED = 299792458;       % WGS-84 value in m / s
EARTH_RATE = 7.2921151467e-5;  % WGS-84 value in rad/s 
MU_EARTH = 3.986005e14;        % WGS-84 value in m^3/s^2

% GPS constants
L1_FREQ = 1575.42e6;                      % Hz (1575.42 MHz)
L1_WAVELENGTH = LIGHT_SPEED / L1_FREQ;

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_pr=[]; prn=[]; pr=[]; pr_orb=[]; pr_errors=[]; obscure_info=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(6,12,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PSEUDO_R for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

if nargin < 7
  model = [1 1 1 1 1 1 1 0 0 0 0 0];     % set to the default value
end % if nargin > 5

if nargin < 8
  seed = 0;     % set to the default value
end % if nargin > 6

if nargin < 9  
  code_noise = 1;
  carrier_noise = .01;
end % if nargin < 8

if nargin < 10  
  carrier_noise = .01;
end % if nargin < 9

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'PSEUDO_R';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_user';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_user;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'x_user';
  estruct.variable(2).req_dim = [901 3; 901 4; 901 5; 1 3; 1 4; 1 5];
  estruct.variable(2).var = x_user;
  
  estruct.variable(3).name = 'v_user';
  estruct.variable(3).req_dim = [901 3; 1 3];
  estruct.variable(3).var = v_user;
  
  estruct.variable(4).name = 't_gps';
  estruct.variable(4).req_dim = [902 2];
  estruct.variable(4).var = t_gps;
  estruct.variable(4).type = 'GPS_TIME';
  
  estruct.variable(5).name = 'x_gps';
  estruct.variable(5).req_dim = [902 4];
  estruct.variable(5).var = x_gps;
  
  estruct.variable(6).name = 'v_gps';
  estruct.variable(6).req_dim = [902 3];
  estruct.variable(6).var = v_gps;
  
  estruct.variable(7).name = 'model';
  estruct.variable(7).req_dim = [1 11];
  estruct.variable(7).var = model;
  
  estruct.variable(8).name = 'seed';
  estruct.variable(8).req_dim = [1 1];
  estruct.variable(8).var = seed;
  
  estruct.variable(9).name = 'code_noise';
  estruct.variable(9).req_dim = [1 1];
  estruct.variable(9).var = code_noise;
  
  estruct.variable(10).name = 'carrier_noise';
  estruct.variable(10).req_dim = [1 1];
  estruct.variable(10).var = carrier_noise;
  
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

clear estruct
%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% Determine how many receivers and how many antenna per receiver are being used
insert_comb_nums = 0;     % no combination numbers
if size(x_user,2) == 4              % more than one receiver is in use
  rec_nums = x_user(:,1); 
  ant_nums = ones(size(x_user,1),1);

% If there are numtiple antenna in use, store the antenna numbers in a separate
% matrix and remove that information from the x_user matrix.  Also we must 
% compute a unique ID number for each receiver/antenna for use in LOS.
% Initialize the flag indicating if combination receiver/ant numbers are used,
elseif size(x_user,2) == 5
  rec_nums = x_user(:,1); 
  ant_nums = x_user(:,2);                                             
  
  % Compute unique receiver/antenna numbers for LOS.  This is limited to 10000 
  % receivers at a time.  If you need to do more than that increase this
  % to something larger and call technical support to let them know.  Please.
  comb_num = x_user(:,1) * 10000 + ant_nums;
  insert_comb_nums = 1;     % combination numbers
    
  x_user = [comb_num x_user(:,3:5)];
else      % input is nx3
  rec_nums = ones(size(x_user,1),1);
  ant_nums = ones(size(x_user,1),1);
  x_user = [rec_nums x_user];      % this sets x_user to nx4 from here on
end % if size(x_user,2) > 3
  
% Compute line-of-sight vectors.  With the vectorized LOS routine, it will
% accept mutiple receiver numbers without an external loop.  If the x_user and
% x_gps matrices are the same length, they are assumed to already be aligned
% and only those PR measurements are computed.
if size(x_user,1) == size(x_gps,1)
  los = x_gps(:,2:4) - x_user(:,2:4);
  t_pr = t_gps;
  pr_ndex = [1:size(x_user,1); 1:size(x_user,1)]';
else
  [t_pr,los,pr_ndex,obscure_info] = los(t_user,x_user,t_gps,x_gps);  
end % if size(x_user,1) == size(x_gps,1)

% Compute the magnitude and direction of the LOS vector.  The magnitude is
% the starting point for the pseudorange (pr) measurement.
[los_norm, pr] = normvect(los); 

% Model Earth rotation
if model(10) == 1
  % Compute the pseudorange in seconds
  pr_seconds = pr ./ LIGHT_SPEED;

  % Correct the user position for Earth rotation, per ICD-GPS-200C-002
  % 25 Sept, 1997, section 20.3.3.4.3.4, page 102.
  earth_rot_corr(:,1) = -EARTH_RATE * x_user(pr_ndex(:,1),2) .* pr_seconds;
  earth_rot_corr(:,2) = EARTH_RATE * x_user(pr_ndex(:,1),1) .* pr_seconds;

  % Use the indices from LOS to obtain which x_user are associated with the LOS,
  % and therefore, PR measurement.
  x_new(:,1:2) = x_user(pr_ndex(:,1),2:3) + earth_rot_corr;
  x_new(:,3) = x_user(pr_ndex(:,1),4);

  % Recompute the LOS and PR with the new user positions
  los = x_gps(pr_ndex(:,2),2:4) - x_new;
  [los_norm, pr] = normvect(los); 

end % if model(10) == 1

% Model satellite motion.  Reference: 'GPS Theory and Practice", 
% Hoffman-Wellenhoff, pages 89-90.
if model(8) == 1
  % Compute the pseudorange in seconds
  pr_seconds = pr / LIGHT_SPEED;

  % compute the relative velocity of the user trajectory and the GPS satellite
  vel = v_gps(pr_ndex(:,2),1:3) - v_user(pr_ndex(:,1),1:3);

  % compute the Doppler shift as the dot product of the LOS with the velocity
  doppler = dot(los_norm',vel')';

  % Compute the correction for satellite motion
  sat_motion = doppler .* pr_seconds;
  
  % Correct the PR measurements
  pr = pr + sat_motion;

else
  sat_motion = zeros(size(pr));  
end % if model(8) == 1

% If a combination was used, replace the combination number with the receiver
% number and regain the original input (less the antenna number)
if insert_comb_nums == 1
  x_user = [rec_nums x_user(:,2:4)]
end % if insert_comb_nums == 1

% rename the GPS satellite number (prn), receiver number, and antenna number
% variable for easier reading later in the code.  Do the same thing for the
% x_user data so it can be easily used in the tropo and iono models where
% azimuth and elevation relative to local level are required.
prn = x_gps(pr_ndex(:,2),1); 
ant_nums = ant_nums(pr_ndex(:,1));
x_user = x_user(pr_ndex(:,1),:);
rec_nums = x_user(:,1);

% Initialize the SA epsilon and dither error matrices and the receiver clocks
sa_eps_err = zeros(size(pr,1),1);
sa_clk_err = zeros(size(pr,1),1);
sa_clk_drift = zeros(size(pr,1),1);
clk_bias = zeros(size(pr,1),1);
clk_drift = zeros(size(pr,1),1);
line_bias = zeros(size(pr,1),1);

% Compute the range portion of the accumulated carrier phase measurement.
% Keep the phase in meters until all computations are complete, then
% convert to phase using the L1 wavelength.
cph = pr;

% Loop over the receivers to generate S/A and receiver clock errors
rec_num_unique = unique(rec_nums);
num_receivers = length(rec_num_unique);

for rec_loops = 1:num_receivers
  % Find all of the receivers with this receiver number
  I = find(x_user(:,1) == rec_num_unique(rec_loops));
  
  % Loop over the number of antenna for this receiver and add the line biases
  % See Parkinson, Vol. 2, page 28, eq 11.  The line bias is applied to each
  % antenna and would appear as a clock bias for that antenna.  However, when
  % multiple antenna are in use, the line bias is not common between the two 
  % antenna unless the cable length is the same.  For this case, set the lin
  % bias sigma to 0 to apply no line bias.
  ant_unique = unique(ant_nums(I));
  for ant_n = 1:length(ant_unique)
    I_this_ant = find(ant_nums == ant_unique(ant_n) & ...
                      x_user(:,1) == rec_num_unique(rec_loops));
    
    % Add the line bias if modeled
    if model(7) == 1 
      % Seed the random number generator such that is is repeatable,
      % but such that all of the line biases are not the same.
      if num_receivers > 1
        if seed == 0
          seed = 1;
        end % if num_receivers > 1
        seed = seed * rec_num_unique(rec_loops) * 1000 + ant_unique(ant_n);
      end % if num_receivers > 1                
      
      rand('seed',seed);                               
      
      % Get the line bias for this receiver and antenna
      if exist('model_data')
        if isfield(model_data,'line_bias')
          % Make sure that this is a nx3 using the error checking function
          estruct.variable(1).name = 'model_data.line_bias';
          estruct.variable(1).req_dim = [901 3];
          estruct.variable(1).var = model_data.line_bias;

          % Call the error checking function
          stop_flag = err_chk(estruct);
  
          if stop_flag == 1           
            fprintf('Invalid line bias inputs to %s.  Using default = 1 meter.\n\n', ...
                     estruct.func_name);
            line_bias_sigma = 1;
          else
         
            I = find(model_data.line_bias(:,1) == rec_num_unique(rec_loops) & ...
                     model_data.line_bias(:,2) == ant_unique(ant_n));
                 
            if ~isempty(I)
              line_bias_sigma = model_data.line_bias(I(1),3);
            else
              line_bias_sigma = 1;
            end % if ~isempty(I) 
          
          end % if stop_flag == 1
        
        else % if isfield(model_data,'line_bias')
          line_bias_sigma = 1;
      
        end % if isfield(model_data,'line_bias')
      
      else
        line_bias_sigma = 1;
      end % if exist('model_data')
        
      line_bias(I_this_ant) = ones(length(I_this_ant),1) * rand(1,1) * line_bias_sigma;  
      pr(I_this_ant) = pr(I_this_ant) + line_bias(I_this_ant);
      cph(I_this_ant) = cph(I_this_ant) + line_bias(I_this_ant);
    else
      line_bias(I_this_ant) = zeros(size(I_this_ant));
    end % if model(7) == 1  
          
  end % for ant_n = 1:length(ant_unique)
  
  % add sa epsilon effect
  if model(1) == 1       
    if exist('model_data')
      if isfield(model_data,'sa_eps')
        % Make sure that this is a 1x1 using the error checking function
        estruct.variable(1).name = 'model_data.sa_eps';
        estruct.variable(1).req_dim = [1 1];
        estruct.variable(1).var = model_data.sa_eps;

        % Call the error checking function
        stop_flag = err_chk(estruct);
  
        if stop_flag == 1           
          fprintf('Invalid SA epsilon inputs to %s.  Using default = 1 meter.\n\n', ...
                   estruct.func_name);
          % Set the sigma value for the SA epsilon model to the RTCM default of 23
          sa_eps_sigma = 23;    % meters                                          
        else
          % Set the sigma value for the SA epsilon model to the RTCM default of 23
          sa_eps_sigma = 23;    % meters                                          
        end % if stop_flag == 1
    
      else
        % Set the sigma value for the SA epsilon model to the RTCM default of 23
        sa_eps_sigma = 23;    % meters                                          
      end % if isfield(model_data,'sa_eps')
    
    else
      % Set the sigma value for the SA epsilon model to the RTCM default of 23
      sa_eps_sigma = 23;    % meters                                          
    end % if exist('model_data')

    % compute the SA sa_eps errors and apply them to the PRs 
    % all done in the sa_eps function
    [pr_1, sa_eps_err_1] = sa_eps([prn(I) pr(I)], sa_eps_sigma, seed);
    
    pr(I) = pr_1;
    sa_eps_err(I) = sa_eps_err_1;

    % apply the SA epsilon error to the phase data also
    cph(I) = cph(I) + sa_eps_err(I);
  
  % Insert user supplied epsilon model here
  elseif model(1) == 2
    dbstack
    fprintf('Warning from PSEUDO_R.\n');
    fprintf('No user supplied epsilon supplied.  \n');
    fprintf('Setting epsilon contribution to PR to zero (0).\n');
    sa_eps_err(I) = zeros(size(pr(I),1),1);  

  else
    sa_eps_err(I) = zeros(size(pr(I),1),1);  
  end % model(1) == 1  

  % add dither effect
  if model(2) == 1 
    if exist('model_data')
      if isfield(model_data,'sa_clock')
        % Make sure that this is a 1x1 using the error checking function
        estruct.variable(1).name = 'model_data.sa_clock';
        estruct.variable(1).req_dim = [1 3; 3 1];
        estruct.variable(1).var = model_data.sa_clock;

        % Call the error checking function
        stop_flag = err_chk(estruct);
  
        if stop_flag == 1           
          fprintf('Invalid SA dither (clock) inputs to %s.  Using default model.\n\n', ...
                   estruct.func_name);
          % set the sigma value for the SA clock model to the RTCM default values
          sa_model_data = [23 .28 118];    % sigma PR, sigma PR-rate, tau 
                                         % (decorrelation time) (m, m/s, dimensionless)                                          
        else
          % set the sigma value for the SA clock model to the RTCM default values
          sa_model_data = [23 .28 118];    % sigma PR, sigma PR-rate, tau 
                                         % (decorrelation time) (m, m/s, dimensionless)                                          
        end % if stop_flag == 1
    
      else
        % set the sigma value for the SA clock model to the RTCM default values
        sa_model_data = [23 .28 118];    % sigma PR, sigma PR-rate, tau 
                                       % (decorrelation time) (m, m/s, dimensionless)                                          
      end % if isfield(model_data,'sa_clock')
    
    else
      % set the sigma value for the SA clock model to the RTCM default values
      sa_model_data = [23 .28 118];    % sigma PR, sigma PR-rate, tau 
                                       % (decorrelation time) (m, m/s, dimensionless)                                          
    end % if exist('model_data')
    
    % compute the dither effects using the default model (RTCA parameters)
    [sa_clk_err(I), sa_clk_drift(I)] = sa_clock(t_pr(I,:), prn(I), sa_model_data, seed);

    pr(I) = pr(I) + sa_clk_err(I); 
    cph(I) = cph(I) + sa_clk_err(I);                 
  
  % Insert user supplied SA model here
  elseif model(2) == 2
    dbstack
    fprintf('Warning from PSEUDO_R.\n');
    fprintf('No user supplied SA model supplied.  \n');
    fprintf('Setting SA contribution to PR to zero (0).\n');
    sa_clk_err(I) = zeros(size(pr(I),1),1);  

  else
    sa_clk_err(I) = zeros(size(pr(I),1),1);  
  end % model(2) == 1  

  % Add random ambiguities to the phase data (parameter N in Parkinson, page 8,
  % equation 4).
  [cph(I), ambig(I)] = addambig([prn(I) cph(I)]);

  % compute user clock bias
  if model(5) == 1
    if exist('model_data')
      if isfield(model_data,'rec_clock')
        % Make sure that this is a 1x2 using the error checking function
        estruct.variable(1).name = 'model_data.rec_clock';
        estruct.variable(1).req_dim = [1 2];
        estruct.variable(1).var = model_data.rec_clock;

        % Call the error checking function
        stop_flag = err_chk(estruct);
  
        if stop_flag == 1           
          fprintf('Invalid receiver clock inputs to %s.  Using default model.\n\n', ...
                   estruct.func_name);
          % Set the value for the clock model to the default
          clk_model = [4e-19 1.58e-18];     % bias and frequency noise values
        else
          % Set the value for the clock model to the default
          clk_model = [4e-19 1.58e-18];     % bias and frequency noise values
        end % if stop_flag == 1
    
      else
        % Set the value for the clock model to the default
        clk_model = [4e-19 1.58e-18];     % bias and frequency noise values
      end % if isfield(model_data,'rec_clock')

    else
      % Set the value for the clock model to the default
      clk_model = [4e-19 1.58e-18];     % bias and frequency noise values
    end % exist('model_data')
    
    % Compute the clock bias using the model. Use the receiver number in
    % conjunction with the provided seed to generate a new seed value.  
    % This will keep all of the clocks from the different receivers from
    % being the same.    
    if num_receivers > 1
    if seed == 0
      seed = 1;
    end % if num_receivers > 1
    seed = seed * rec_num_unique(rec_loops);
    end % if num_receivers > 1
    
    [clk_bias(I) clk_drift(I)] = clockerr(t_pr(I,:), clk_model, seed);

    % add the clock bias to the pseudo-range and accumulated phase measurements
    pr(I) = pr(I) + clk_bias;
    cph(I) = cph(I) + clk_bias;

  % Insert user supplied receiver clock model here
  elseif model(5) == 2
    dbstack
    fprintf('Warning from PSEUDO_R.\n');
    fprintf('No user supplied receiver clock model supplied.  \n');
    fprintf('Setting receiver clock contribution to PR to zero (0).\n');
    clk_bias(I) = zeros(size(pr(I),1),1);  
    clk_drift(I) = zeros(size(pr(I),1),1);  
    
  else
    clk_bias(I) = zeros(size(pr(I),1),1);  
    clk_drift(I) = zeros(size(pr(I),1),1);  
  end % if model(5) == 1

end % for rec_loops = 1:num_receivers

% Begin the processing for the atmospheric models (troposphere and ionosphere).
% Compute azimuth and elevation in the NED frame.  This is how the tropo and iono
% models want the data.
ref_lla = ecef2lla(x_user(:,2:4));
[x_gps_ned] = ecef2ned(los, ref_lla);
[az, el] = ned2azel(x_gps_ned);     
  
% add troposphere
if model(3) == 1
  
  % Check if there is troposphere model data supplied
  if exist('model_data')
    if isfield(model_data,'tropo')
      % Make sure that this is a 1x3 using the error checking function
      estruct.variable(1).name = 'model_data.tropo';
      estruct.variable(1).req_dim = [1 3];
      estruct.variable(1).var = model_data.tropo;

      % Call the error checking function
      stop_flag = err_chk(estruct);
  
      if stop_flag == 1           
        fprintf('Invalid tropo model inputs to %s.  Using default model.\n\n', ...
                 estruct.func_name);
        % Set the value for the tropo model to the default
        trop_model = [1013.25 288.15 11.691];    
      else
        % Set the value for the tropo model to the default
        trop_model = [1013.25 288.15 11.691];    
      end % if stop_flag == 1
    
    else
      % Set the value for the tropo model to the default
      trop_model = [1013.25 288.15 11.691];    
    end % if isfield(model_data,'tropo')
  
  else
    % Set the value for the tropo model to the default
    trop_model = [1013.25 288.15 11.691];    
  end % exist('model_data')
 
  % Find all of the elevation angles above to to compute iono for
  trop_dry = zeros(size(pr));
  trop_wet = zeros(size(pr));
  
  I = find(el > 0);
  
  if ~isempty(I)
    % Compute the dry and wet troposphere components
    [trop_dry_mdl, trop_wet_mdl] = tropdlay(el(I),trop_model);           
    trop_dry(I) = trop_dry_mdl;
    trop_wet(I) = trop_wet_mdl;
  end % if ~isepmty(I)
    
  % add the troposphere error to the pr
  pr = pr + trop_dry + trop_wet;

  % add the troposphere error to the cph
  cph = cph + trop_dry + trop_wet;
  
% Insert user supplied troposphere model here
elseif model(3) == 2
  dbstack
  fprintf('Warning from PSEUDO_R.\n');
  fprintf('No user supplied troposphere model supplied.  \n');
  fprintf('Setting ionosphere contribution to PR to zero (0).\n');
  trop_wet = zeros(size(az,1),1);
  trop_dry = zeros(size(az,1),1);

else
  trop_wet = zeros(size(az,1),1);
  trop_dry = zeros(size(az,1),1);

end % if model(3) == 1  

% add ionosphere
if model(4) == 1            
  % Check if there is troposphere model data supplied
  if exist('model_data')
    if isfield(model_data,'iono')
      % Make sure that this is a 1x8 using the error checking function
      estruct.variable(1).name = 'model_data.iono';
      estruct.variable(1).req_dim = [1 8];
      estruct.variable(1).var = model_data.iono;

      % Call the error checking function
      stop_flag = err_chk(estruct);
  
      if stop_flag == 1           
        fprintf('Invalid iono model inputs to %s.  Using default model.\n\n', ...
                 estruct.func_name);
        % Set the value for the tropo model to the default
        iono_params = [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
                       0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06]; 
      else
        % Set the value for the tropo model to the default
        iono_params = [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
                       0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06]; 
      end % if stop_flag == 1
    
    else
      % Set the value for the tropo model to the default
      iono_params = [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
                     0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06]; 
    end % if isfield(model_data,'iono')
  else
    % Set the value for the tropo model to the default
    iono_params = [0.1211D-07  0.1490D-07 -0.5960D-07 -0.1192D-06 ...
                   0.9626D+05  0.8192D+05 -0.1966D+06 -0.3932D+06]; 
  end % if exist('model_data')
  
  % Find all of the elevation angles above to to compute iono for
  iono = zeros(size(pr));
  
  I = find(el > 0);
  if ~isempty(I)
    iono_dlay = ionodlay(t_pr(I,:),ref_lla(I,:),az(I),el(I),iono_params);
    iono(I) = iono_dlay;
  end % if ~isepmty(I)
  
  % add the ionosphere error to the pr
  pr = pr + iono;
  
  % add the ionosphere error to the cph (the sign is opposite from
  % the pr because the ionosphere delays the code (pr) and advances
  % the carrier (cph))
  cph = cph - iono;

% Insert user supplied ionosphere model here
elseif model(4) == 2                  
  dbstack
  fprintf('Warning from PSEUDO_R.\n');
  fprintf('No user supplied ionosphere model supplied.  \n');
  fprintf('Setting ionosphere contribution to PR to zero (0).\n');
  iono = zeros(size(az,1),1);

else
  iono = zeros(size(az,1),1);

end % if model(4) == 1

% Add GPS satellite clock model from ephemeris
if model(9) == 1            
  % Check if there is ephemeris data supplied
  if exist('ephem')
    % Make sure that this is a nx24 using the error checking function
    estruct.variable(1).name = 'ephem';
    estruct.variable(1).req_dim = [901 24];
    estruct.variable(1).var = ephem;

    % Call the error checking function
    stop_flag = err_chk(estruct);
  
    if stop_flag == 1           
      fprintf('Invalid ephemeris inputs to %s.  ',estruct.func_name);
      fprintf('No GPS satellite clock correction applied.\n\n');
      % Set the value for the tropo model to the default
      sat_clock = zeros(size(az,1),1);
    else
      % Compute the satellite clock correction for each GPS satellite
      % based on the ephemeris data.
      % Start by find the unique satellite IDs and initialing the sat_clock
      sat_ids = unique(ephem(:,1));
      sat_clock = zeros(size(az,1),1);
      
      % Compute the time past the Toc (time of clock) parameter.
      pr_time = gpst2sec(t_pr);
      
      % Loop over the valid satellite IDs and compute the clock correction
      % for that satellite.
      for i = 1:length(sat_ids)
        % Find all of the PR and ephem with this satellite ID
        I_pr = find(prn == sat_ids(i));
        I_ephem = find(ephem(:,1) == sat_ids(i));
        
        % Find the time past the Toc (time of clock)
        toc = gpst2sec([ephem(I_ephem(1),19) ephem(I_ephem(1),20)]);
        
        % Compute the time between toc and the PR times
        delta_t = pr_time(I_pr) - toc;
        
        % Compute the satellite clock correction to the PR in meters
        sat_clock(I_pr) = (ephem(I_ephem(1),21) + ...
                           ephem(I_ephem(1),22) * delta_t + ...
                           ephem(I_ephem(1),23) * delta_t.^2) * LIGHT_SPEED;
      end % for i = 1:length(sat_ids)
    end % if stop_flag == 1
    
  else
    sat_clock = zeros(size(az,1),1);
  end % if exist('ephem')

  % add the satellite clock error to the pr
  pr = pr + sat_clock;
  
  % add the ionosphere error to the cph (the sign is opposite from
  % the pr because the ionosphere delays the code (pr) and advances
  % the carrier (cph))
  cph = cph + sat_clock;

else
  sat_clock = zeros(size(az,1),1);

end % if model(9) == 1

% Add relativity correction from ephemeris data per IC-GPS-200C, 10 Oct 1993
% pages 88-89 section 20.3.3.3.3.1, User Algorithm for the SV Clock Correction.
if model(11) == 1            
  % Check if there is ephemeris data supplied
  if exist('ephem')
    % Make sure that this is a nx24 using the error checking function
    estruct.variable(1).name = 'ephem';
    estruct.variable(1).req_dim = [901 24];
    estruct.variable(1).var = ephem;

    % Call the error checking function
    stop_flag = err_chk(estruct);
  
    if stop_flag == 1           
      fprintf('Invalid ephemeris inputs to %s.  ',estruct.func_name);
      fprintf('No relativistic corrections applied.\n\n');
      % Set the value for the relativistic effects to zero
      relative = zeros(size(az,1),1);
    else                     
      % Define the constant F, per ICD-GPS-200, section 20.3.3.3.3.1
      F = -4.442807633e-10;
      
      % Compute the PR times in seconds past GPS epoch
      pr_time = gpst2sec(t_pr);

      % Start by find the unique satellite IDs and initialing the sat_clock
      sat_ids = unique(ephem(:,1));
      relative = zeros(size(az,1),1);
      
      % Loop over the valid satellite IDs and compute the relativistic 
      % correction for that satellite.
      for i = 1:length(sat_ids)
        % Find all of the PR and ephem with this satellite ID
        I_pr = find(prn == sat_ids(i));
        I_ephem = find(ephem(:,1) == sat_ids(i));
        
        % Find the Toe (time of ephemeris) in GPS seconds
        toe = gpst2sec([ephem(I_ephem(1),19) ephem(I_ephem(1),17)]);
        
        % Compute the time between toc and the PR times
        delta_t = pr_time(I_pr) - toc;
         
        % Compute some orbit terms
        a12 = ephem(I_ephem(1),5);  % get sqrt(a) = ephem(:,5) (m)
        a = ephem(I_ephem(1),5).^2; % get a = ephem(:,5)^2 (m)
        e = ephem(I_ephem(1),4);    % load eccentricity to a local variable
        n0 = sqrt(MU_EARTH / a.^3);    % compute the nominal mean motion (rad/s)
        n = n0 + ephem(I_ephem(1),3);     % corrected mean montion
                                              % delta_n = ephem(:,3) 
        % Compute the mean anomaly (M)
        M = ephem(I_ephem(1),2) + n .* delta_t; % M0 = ephem(:,2)   
        
        % Compute the eccentric anomaly
        [E] = keplr_eq(M,e);
        
        % Compute the relativistic correction to the PR in meters
        relative(I_pr) = (F .* e .* a12 .* sin(E)) * LIGHT_SPEED;
        
      end % for i = 1:length(sat_ids)
    end % if stop_flag == 1
    
  else
    % Set the value for the relativistic effects to zero
    relative = zeros(size(az,1),1);
  end % if exist('ephem')
    
  % add the satellite clock error to the pr
  pr = pr + relative;
  
  % add the ionosphere error to the cph (the sign is opposite from
  % the pr because the ionosphere delays the code (pr) and advances
  % the carrier (cph))
  cph = cph + relative;

else
  relative = zeros(size(az,1),1);

end % if model(11) == 1

% add receiver white noise 
if model(6) == 1
  % Compute the clock bias using the model. Use the receiver number in
  % conjunction with the provided seed to generate a new seed value.  
  % This will keep all of the clocks from the different receivers from
  % being the same.    
  if num_receivers > 1
    if seed == 0
      seed = 1;
    end % if num_receivers > 1
    seed = seed * rec_num_unique(rec_loops);
  end % if num_receivers > 1

  % set the seed value of the random number generator 
  randn('seed',seed);
  
  % compute white noise for the code and carrier
  code_white = randn(size(pr,1),1) * code_noise;
  carrier_white = randn(size(pr,1),1) * carrier_noise;

  pr = pr + code_white; 
  cph = cph + carrier_white;
   
else           

  code_white = zeros(size(pr,1),1);
  carrier_white = zeros(size(pr,1),1);
end % if model(6) == 1

% convert the phase measurements from meters to phase
cph = cph / L1_WAVELENGTH;

% build up the PR error matrix for output
pr_errors = [sa_eps_err sa_clk_err trop_wet trop_dry iono clk_bias clk_drift ...
             code_white carrier_white line_bias sat_motion sat_clock relative];

% add the phase data to the output variable pr
pr = [pr cph];

%%%%% END ALGORITHM CODE %%%%%

% end of PSEUDO_R
