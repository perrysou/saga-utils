function [t_dpr, dpr] = diffcorr(t_pr, pr, pr_orb, pr_orb_vel, base, base_model, ephem);

% [t_dpr, dpr] = diffcorr(t_pr, pr, pr_orb, pr_orb_vel, base, base_model, ephem);
%
% Function to compute differential pseudorange (PR) corrections (DPR) 
% at a base station with known coordinates.
%
% Input:
%   t_pr    - GPS time associated with the PR collected at a base station
%              [GPS_week GPS_sec] (nx2)
%              valid GPS_week values are 1-3640 (years 1980-2050)
%              valid GPS_sec values are 0-604799
%   pr     - satellite number and PR for base station pseudo-ranges (nx2)
%             [prn pr] PR is in meters or (nx3) if using carrier phase
%             data [prn pr cph] where chp is in cycles
%   pr_orb  - GPS/GLONASS satellite positions as a function of time
%             associated with this base station PR (nx3) (meters)   
%   pr_orb_vel - GPS/GLONASS satellite velocities at t_pr times (nx3) (m/s)
%   base    - base station vector location (1x3) (meters)
%   base_model  - flags controlling which contributions to the 
%                  PR errors are modeled (optional).  (1x11)
%                [sa_eps dither troposphere ionosphere receiver_clock 
%                 receiver_noise line_bias sat_motion sat_clock earth_rotation
%                 relativity]
%                a value of 1 indicates useage of the model
%                and a value of zero indicates no use of that model.  Use values
%                of 2 to implement user supplied models.  See the code for where
%                to insert the user models.  A warning is given if a user model
%                is selected and none is supplied.
%                Default = [0 0 0 0 0 0 0 0 0 0 0]. 
%   ephem       - ephemeris matrix for all satellites (nx24). (optional)
%                  Used to compute satellite clock. If not provided, no 
%                  GPS satellite clock effects will be computed
%                  The columns of ephemeris matrix are ...
%                  [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                  i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                  GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                  Ephemeris parameters are from ICD-GPS-200 with the 
%                  exception of perigee_rate. 
%
%   Note: Coordinate systems for the base and pr_orb data must be the same.
%         Either both in ECEF or ECI frame in meters.
% Output:
%   t_dpr   - GPS time associated with the differential correction, 
%              [GPS_week GPS_sec] (nx2)
%   dpr     - satellite number and differential correction (nx2)
%              [prn_diff dpr] correction (dpr) in meters
%
% See also ADD_DPR, PSEUDO_R, SA_CLOCK, CLOCKERR 

% Written by: Jimmy LaMance 4/3/97 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, LSNAV, PSEUDO_R

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_dpr=[]; dpr=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(5,7,nargin);
if ~isempty(msg)
  fprintf('%s  See help on DIFFCORR for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

if nargin < 6
  base_model = [0 0 0 0 0 0 0 0 0 0 0];
end % if nargin < 6

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'DIFFCORR';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_pr;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'pr';
  estruct.variable(2).req_dim = [901 2; 901 3];
  estruct.variable(2).var = pr;
  
  estruct.variable(3).name = 'pr_orb';
  estruct.variable(3).req_dim = [901 3];
  estruct.variable(3).var = pr_orb;
  
  estruct.variable(4).name = 'base';
  estruct.variable(4).req_dim = [1 3];
  estruct.variable(4).var = base;
  
  estruct.variable(5).name = 'base_model';
  estruct.variable(5).req_dim = [1 11];
  estruct.variable(5).var = base_model;

  estruct.variable(6).name = 'pr_orb_vel';
  estruct.variable(6).req_dim = [901 3];
  estruct.variable(6).var = pr_orb_vel;

  if nargin == 7
    estruct.variable(7).name = 'ephem';
    estruct.variable(7).req_dim = [902 24];
    estruct.variable(7).var = ephem;
  end % if nargin == 7
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% allocate output variable
t_dpr = ones(size(t_pr)) * inf;
dpr = ones(size(pr)) * inf;

% solve for base station clock bias
% this is a crude clock bias solution, but it doesn't matter for DGPS 
% because the base station clock solution is not separable from the
% remote clock solution.  this is used to minimize the values for each
% of the corrections
[t_nav,x_nav,num_sats,nav_index] = lsnav(t_pr,pr(:,2),[pr(:,1) pr_orb],[base 0]);

% find the solutions with 4 or more satellite
I4 = find(num_sats >= 4);
Iless4 = find(num_sats < 4);

% compute a clock bias with the more than 4 solutionss
clock_bias(Iless4) = zeros(size(Iless4,1),1);
clock_bias(I4) = x_nav(I4,4);

% expand the base station location to be nx3 (instead of 1x3)
base = ones(size(t_pr,1),1) * base;

% Compute the 'true' pr and cph measurements that should have been observed
% at the base station using the pseudo_r function with no error modeling.
% The RTCM standard reccommends that no error modeling be applied to base
% station corrections.
seed = 0;                  % random number generator seed 
code_noise = 0;            % don't apply any noise to the code measurement
carrier_noise = 0;         % don't apply any carrier noise to the measurements
v_base = zeros(size(base));

if nargin == 7
  [t_pr_ref,prn_ref,pr_ref,pr_errors_ref] = ...
     pseudo_r(t_pr, base, v_base, t_pr, [pr(:,1) pr_orb], pr_orb_vel,...
              base_model, seed, code_noise, carrier_noise,ephem);
else
  [t_pr_ref,prn_ref,pr_ref,pr_errors_ref] = ...
     pseudo_r(t_pr, base, v_base, t_pr, [pr(:,1) pr_orb], pr_orb_vel,...
              base_model, seed, code_noise, carrier_noise);
end % if nargin == 7

% fill in the satellite number portion of the correction matrix
dpr(:,1) = prn_ref;                       % satellite numbers

% if working with carrier phase data
if size(dpr,2) == 3
  dpr(:,2:3) = pr_ref - pr(:,2:3) + clock_bias(nav_index)' * ones(1,2);
  
  % remove 'most' of the ambiguities from the carrier phase corrections
  % using the initial carrier phase ambiguities        
  % start by sorting through and figuring out which satellites are visibile
  prn_sort = sort(dpr(:,1));
  prn_change = find(diff(prn_sort) ~= 0);

  % create a matrix that has sorted and reduced prns [1 2 4 5 6 8 ... 28]
  prn = [prn_sort(prn_change); prn_sort(length(prn_sort))]; 
  
  % compute the total number of prns 
  num_prns = length(prn);

  % loop over the prns that are visible and remove the ambiguities based on the
  % first available data from that prn
  for i = 1:num_prns
    I = find(dpr(:,1) == prn(i));
    if any(I)
      dpr(I,3) = dpr(I,3) - fix(dpr(I(1),3));
    end % if any(I)
  end % for i = 1:num_prns
  
else   % there is only PR measurements passed into the diffcorr function
  dpr(:,2) = pr_ref(:,1) - pr(:,2) + clock_bias';

end % if size(dpr,2) == 3
  
% assign output variables
t_dpr = t_pr_ref;                         % time output

%%%%% END ALGORITHM CODE %%%%%

% end of DIFFCORR
