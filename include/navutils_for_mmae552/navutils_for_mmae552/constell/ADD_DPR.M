function [t_cpr,prn_cpr,cpr,cpr_index] = ...
                               add_dpr(t_pr,pr,t_dpr,dpr,dprr,max_latency);

% [t_cpr,prn_cpr,cpr,cpr_index] = add_dpr(t_pr,pr,t_dpr,dpr,dprr,max_latency);
%
% Function to apply differential corrections to pseudo-range (PR) and accumulated
% carrier phase (CPH) measurements.
%
% Input:
%   t_pr  - time associated with the PR, [GPS_week GPS_sec] (nx2)
%            valid GPS_week values are 1-3640 (years 1980-2050)
%            valid GPS_sec values are 0-604799
%   pr    - satellite number and PR (and CPH) for each t_pr [prn pr] (nx2) or 
%            [prn pr cph] (nx3) pr is in meters and cph is in cycles
%            (or same units as dpr)
%   t_dpr - time associated with the differential correction (DPR), 
%            [GPS_week GPS_sec] (kx2)
%   dpr   - satellite number and DPR (and dCPH) for each t_dpr [dprn dpr] (kx2) 
%            or [dprn dpr dcph] (kx3)  dpr units are meters and dcph is in cycles 
%            (or same units as pr)
%   dprr  - differential correction rate associated with the 
%            corresponding dpr and t_dpr (kx1) (optional) 
%            default: dprr = 0, rate is in meters/s (or same units 
%            as pr per second.  
%   max_latency - maximum latency for differential corrections to be applied (s)
%                  (optional), default = 10
% Output:
%   t_cpr   - time associated with the corrected pseudo-range (CPR), 
%              [GPS_week GPS_sec] (mx2)
%   prn_cpr - satellite number for base station CPRs (mx1)
%   cpr     - pseudo-range correction associated with the 
%              corresponding t_cpr and prn_cpr (mx1) or (mx2) if carrier phase
%              corrections are computed
%   cpr_index - index matrix that relates the corrected pseudo-ranges
%                to the input pr measurements.
%
% See also DIFFCORR, PSEUDO_R, SA_CLOCK, CLOCKERR 

% Written by: Jimmy LaMance 4/7/97 
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
t_cpr=[]; prn_cpr=[]; cpr=[]; cpr_index=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(4,6,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ADD_DPR for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check if a differential range rate correction is provided
if nargin < 5
  dprr = zeros(size(dpr,1),1);
end % if nargin < 6

% check if a maximum latency is provided
if nargin < 6
  max_latency = 10;
end % if nargin < 6

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ADD_DPR';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 't_pr';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_pr;
  estruct.variable(1).type = 'GPS_TIME';
  
  estruct.variable(2).name = 'pr';
  estruct.variable(2).req_dim = [901 2; 901 3];
  estruct.variable(2).var = pr;
  
  estruct.variable(3).name = 't_dpr';
  estruct.variable(3).req_dim = [902 2];
  estruct.variable(3).var = t_dpr;
  
  estruct.variable(4).name = 'dpr';
  estruct.variable(4).req_dim = [902 2; 902 3];
  estruct.variable(4).var = dpr;
  
  estruct.variable(5).name = 'dprr';
  estruct.variable(5).req_dim = [902 1];
  estruct.variable(5).var = dprr;
  
  estruct.variable(6).name = 'max_latency';
  estruct.variable(6).req_dim = [1 1];
  estruct.variable(6).var = max_latency;
  
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

% break out the cph and dcph measurements from the input variables
if size(pr,2) == 3
  cph = pr(:,3);
  dcph = dpr(:,3);
else
  cph = [];               % initialize the cph measurements to blank,
  dcph = [];              % this enables checking on the size later
end % if size(pr,2) == 3

% break out the pr and prn from the input variables,
% the variable names pr and dpr are still used without the prn associated
% with them.  this saves on memory, but decreases the readability for the next
% 4 lines of code.
prn_dpr = dpr(:,1);
dpr = dpr(:,2);

prn_pr = pr(:,1);
pr = pr(:,2); 

% match times and prn numbers for corrections and raw pseudo-ranges 
% start by establishing which PRN numbers we have
prn_all = [prn_pr; prn_dpr];
prn_all_sort = sort(prn_all);
diff_prn = diff(prn_all_sort);
I_prn = find(diff_prn ~= 0);
prn = [prn_all_sort(1); prn_all_sort(I_prn+1)]; 

clear I_prn prn_all prn_all_sort diff_prn

num_prns = size(prn,1);

% allocate the output matrices the same size as the differential corrections
% and fill with inf value.  At the end, only the valid corrections will be 
% filled with numbers (other than inf).  Then a subset of the allocated
% matrix will be returned with the indices indicating valid value.
t_cpr = ones(size(t_pr)) * inf;
prn_cpr = ones(size(prn_pr)) * inf;
cpr = ones(size(pr)) * inf;
cpr_index = [];

% now loop over each satellite
for i = 1:num_prns 
  % sort out the PR measurements for this satellite
  I_pr = find(prn_pr == prn(i));
  I_dpr = find(prn_dpr == prn(i));
  
  if any(I_pr) & any(I_dpr)
    % compute linear time for this sv
    t_pr_temp = t_pr(I_pr,1) * 604800 + t_pr(I_pr,2);

    % compute linear time for this sv differential corrections
    t_dpr_temp = t_dpr(I_dpr,1) * 604800 + t_dpr(I_dpr,2);    
  
    % make sure that the first time for the PR measurement is later than
    % or equal to the time for the first differential correction for this
    % satellite

    I_not_early = find(t_pr_temp >= min(t_dpr_temp));
    I_pr = I_pr(I_not_early); 
  
    t_pr_temp = t_pr_temp(I_not_early);
    
    % set up a temporary matrix with the PR and CPH measurements 
    % for this satellite
    pr_temp = pr(I_pr);
    
    if size(cph,1) >= 1
      cph_temp = cph(I_pr);
      dcph_temp = dcph(I_dpr,:); 
    end % if size(cph,1) >= 1
  
    % sort out the differential corrections for this satellite
    dpr_temp = dpr(I_dpr,:); 
  
    % construct a time matrix with all of the temp times, PR, and corrections
    t_all_temp = [t_dpr_temp; t_pr_temp];
    [t_all_sorted, I_sort] = sort(t_all_temp);

    % make sure that the first element in the corrected matrix is a differential
    % correction.  this insures that only corrections for times past will
    % be used to correct PR measurements.
    I1 = find(I_sort < size(t_dpr_temp,1));  
    start_sort = I1(1);
    stop_sort = size(I_sort,1);

    I_sort = I_sort(start_sort:stop_sort);
    t_all_sorted = t_all_sorted(start_sort:stop_sort);
    
    % construct a pr correction matrix of the same form as the time matrix  
    pr_corr_matrix = [dpr_temp; pr_temp];
    pr_corr_matrix = pr_corr_matrix(I_sort);

    % construct a CPH correction matrix of the same form as the time matrix  
    if size(cph,1) >= 1
      cph_corr_matrix = [dcph_temp; cph_temp];
      cph_corr_matrix = cph_corr_matrix(I_sort);
    end % if size(cph,1) >= 1

    % construct a pr rate correction matrix of the same form as pr_corr_matrix
    prr_corr_matrix = [dprr(I_dpr); pr_temp];
    prr_corr_matrix = prr_corr_matrix(I_sort);  
  
    % compute time differences (this has both PR and DPR)   
    dt = diff(t_all_sorted);

    % find the sorted observations to be corrected (they were at the
    % end of the original matrix, now they are interspersed)  This finds
    % the sorted index with indices greater than the number of corrections
    number_of_dprs = size(dpr_temp,1);
    I_t2 = find(I_sort > number_of_dprs);   % 2 corresponds to PRs

    % now search back to find the correction that is closest in time before
    % the PR measurement (I_t2 index)
    % the correction that corresponds to the PR meas. should be just before it
    % if there are no skips in the data
    I_t1 = I_t2 - 1;                        % 1 corresponds to the corrections

    % verify that there are no skips in the PR data.  a skip would put
    % pseudorange measurements next to each other instead of a every
    % other set of corrections and PR measurements
    I_big = find(I_sort(I_t1) > number_of_dprs);  
    
    if any(I_big)
      finished = 0;
      while finished == 0
        I_t1(I_big) = I_t1(I_big) - 1;
        clear I_big
      
        % check again
        I_big = find(I_sort(I_t1) > number_of_dprs);
      
        if ~any(I_big)
          finished = 1;
        end % if ~any(I_big)
      end % while finished == 0
    end % if any(I_big)
      
    % verify that none of the I_t1 variables have a negative index
    if any(find(I_t1) < 1)
      fprintf('Error sorting PR and differential corrections.\n')
      keyboard
    end % if any(find(I_t1) < 1)

    % find the time differences between the corrections and the PR measurements
    % this is the time since the differential correction
    delta_t = t_all_sorted(I_t2) - t_all_sorted(I_t1);
    
    % see if any of the corrections are 'late'
    I_late = find(delta_t > max_latency);    
    
    % if there are late correction, fill the corresponding elements of the 
    % pr_corr_matrix with inf to flag this as an invalid correction
    if any(I_late)
      pr_corr_matrix(I_t1(I_late)) = ones(size(I_late,1),1) * inf;
      clear I_late
    end

    % compute the corrected PR measurement as  
    % PR + correction + PR_correction_rate * delta_t
  
    % corr PR  =   PR     +      correction      + ...
    %              correction rate       *  delta time
    cpr(I_pr,1) = pr_temp + pr_corr_matrix(I_t1) + ...
                  prr_corr_matrix(I_t1) .* delta_t; 

    % compute the corrected CPH measurement as the orginial 
    % CPH + correction
  
    if size(cph,1) >= 1
      % corr CPH  =   CPH    +    correction    
      cpr(I_pr,2) = cph_temp + cph_corr_matrix(I_t1); 
    end % if size(cph,1) >= 1                    

    % output time matrix
    t_cpr(I_pr,:) = t_pr(I_pr,:);
  
    % output prn matrix
    prn_cpr(I_pr,:) = prn_pr(I_pr,:);
    
    clear I_pr I_dpr I_big pr_corr_matrix prr_corr_matrix I_t1 I_t2 dt delta_t
    clear pr_temp dpr_temp t_pr_temp t_dpr_temp
  end % if any(I_pr)
  
end % for i = 1:num_prns  

% find the valid corrected PR measurements
cpr_index = find(cpr(:,1) ~= inf);

% fill in the output time matrix and corrected PR measurement
t_cpr = t_pr(cpr_index,:);
prn_cpr = prn_pr(cpr_index);
cpr = cpr(cpr_index,:);

%%%%% END ALGORITHM CODE %%%%%

% end of ADD_DPR
