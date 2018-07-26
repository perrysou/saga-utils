function [t_out, num_sats, id_num_sats] = num_vis(t_vis, observer_ids);

% [t_out, num_sats, id_num_sats] = num_vis(t_vis, observer_ids);
%
% Computes the number of satellites visible as a function of time.
%
% Input:
%   t_vis   - time that observed satellites are visible in GPS format (nx2)
%             [GPS_week GPS_sec]
%             GPS_week values may range from 0 to 3640
%             GPS_sec values may range from 0 to 604800
%   observer_ids - Optional identification numbers for observers
%                  (i.e. PRN, ground station #) (nx1). If not provided, 
%                  it is assumed that there only 1 observer, and id_num_sats
%                  will be filled with ones.
% Output:
%   t_out       - output time matrix (kx2) [GPS_week GPS_sec]
%   num_sats    - number of satellites visible at each output time (kx1)
%   id_num_sats - id number of the observer for each time step (kx1)
%

% Written by: Jimmy LaMance 9/9/97
% Modified by: Maria Evans 5/5/98
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK, GPST2SEC, SEC2GPST

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
t_out=[]; num_sats=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on NUM_VIS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'NUM_VIS';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.

  estruct.variable(1).name = 't_vis';
  estruct.variable(1).req_dim = [901 2];
  estruct.variable(1).var = t_vis;
  estruct.variable(1).type = 'GPS_TIME';
  
  if nargin == 2,
    estruct.variable(2).name = 'observer_ids';
    estruct.variable(2).req_dim = [901 1];
    estruct.variable(2).var = observer_ids;
  end;

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

if ~isempty(t_vis),
  % First sort by GPS time
  t_lin = gpst2sec(t_vis);
  [data,i_1] = sort(t_lin);

  % if id numbers are provided, sort the data according to id numbers
  if nargin == 2,
    % Now sort the data by id_number
    new_data = [data observer_ids];
    [data,i_2] = sortrows(new_data,2);
  else,
    data = [data ones(size(data,1),2)];
  end;

  % find changes in t (different sets of observations)
  delta_t = diff(data(:,1));

  % Find breaks in id #1
  delta_id = diff(data(:,2));

  % Determine the time step
  i_nonzero = find(delta_t > 0);
  time_step = min(delta_t(i_nonzero));

  if nargin > 1,
    I_break = find(delta_t > 0 | delta_id > 0);
  else,
    I_break = find(delta_t > 0);
  end;

  % fill in number of visible satellites 
  t_sec_1 = [data(I_break,1)];
  num_sats = [diff(I_break); size(data,1) - I_break(length(I_break));];
  id_num_sats = data(I_break,2);

  % Check for time breaks greater than time_step
  % This indicates time when there are no visible satellites
  I_break = find(delta_t > time_step);
  if any(I_break),
    t_sec_1 = [t_sec_1; data(I_break,1)+time_step; data(I_break+1,1)-time_step];
    num_sats = [num_sats; zeros(2*length(I_break)*1,1);];
    id_num_sats = [id_num_sats; zeros(2*length(I_break)*1,1);];
    % Resort the data for time
    [t_sec_1,i_1] = sort(t_sec_1);
    num_sats = num_sats(i_1);
    id_num_sats = id_num_sats(i_1);
  end;

  t_out = sec2gpst(t_sec_1);
end;

%%%%% END ALGORITHM CODE %%%%%

% end of NUM_VIS
