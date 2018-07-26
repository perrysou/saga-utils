function [pass_numbers, pass_times, pass_summary] = ...
                        passdata(gps_time, pass_dt, id_nums, other_vis_data)

% [pass_numbers, pass_times, pass_summary] = passdata(gps_time, pass_dt, 
%           					id_nums, other_vis_data);
%
% Determines pass numbers for data that has passed the masking
% tests in vis_data. Also provide summary information about each pass. 
%
% Input:
%   gps_time = GPS times that targets are visible to observers 
%              [GPS week, GPS_sec] (nx2)
%   pass_dt  = Minimum number of seconds between consecutive data
%              determining whether this data is part of the previous
%              datas visibility pass. (Optional) (1x1) (seconds). Should
%              be a number greater than your propagation step size. If
%              not provided, all passes will be given the pass number 1.
%   id_nums  = Optional identification numbers of observer and target objects.
%              (i.e. PRN, ground station #) for two object (nx2). Used to 
%              differentiate pass numbers based on which objects were involved 
%              in the observation.
%   other_vis_data = Optional array of other visible data i.e. azimuth, elevation, 
%                    range, etc. (nxm) (any units)
% Output:
%   pass_numbers = Pass number of each data point. (nx1)
%   pass_times   = [Pass number, GPS week, GPS sec at start of pass, 
%                   duration of pass (sec), id_num1, id_num2] (jx6)
%   pass_summary = Summary of other_vis_data for each pass (jxmxk) where
%                  j = number of passes
%                  m = number of input variables in other_vis_data
%                  k = 4 [value at start of pass, value at end of pass, 
%                         minimum during pass, maximum during pass]
%                  Units match input units in other_vis_data
%
% See Also: VIS_DATA, VIS_E, VIS_O

% Written by: Maria Evans, April 1998
% Copyright (c) 1998 by Constell, Inc.

% functions called: GPST2SEC, SEC2GPST

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
pass_numbers=[];
pass_times = [];
pass_summary = [];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,4,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PASSDATA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'PASSDATA';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'gps_time';
  estruct.variable(1).req_dim = [901 2;];
  estruct.variable(1).var = gps_time;
  estruct.variable(1).type = 'GPS_TIME';
  
  if nargin >= 2,
    estruct.variable(2).name = 'pass_dt';
    estruct.variable(2).var = pass_dt;
    estruct.variable(2).req_dim = [1 1];

    if nargin >= 3,
      estruct.variable(3).name = 'id_nums';
      estruct.variable(3).var = id_nums;
      estruct.variable(3).req_dim = [901 2];
      
      if nargin == 4,
        estruct.variable(4).name = 'other_vis_data';
        estruct.variable(4).var = other_vis_data;
        estruct.variable(4).req_dim = [901 size(other_vis_data,2)];
      end;
    end;
  end;
 
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

% if id numbers are provided, sort the data according to id numbers

if nargin >= 3,
  % Create a combined ID number to sort by
  numdigits = length(int2str(ceil(max(id_nums(:,2)))));
  multiplier = sprintf('%d',zeros(1,numdigits));
  multiplier = str2num(sprintf('1%s',multiplier));
  id_combined = id_nums(:,1)*multiplier + id_nums(:,2);
else,
  id_combined = ones(size(gps_time,1),1);
end;

% First sort by GPS time
if nargin == 4,
  new_data = [gpst2sec(gps_time) id_combined id_nums other_vis_data];
else,
  new_data = [gpst2sec(gps_time) id_combined id_nums]; 
end;
[data,i_1] = sortrows(new_data,1);

% Now sort the data by id_number
new_data = data;
[data,i_2] = sortrows(new_data,2);

% if pass_dt is provided
if nargin >= 2,
  % Check for breaks in data
  delta_t = diff(data(:,1));
  delta_id = diff(data(:,2));
  I_break = find(delta_t > pass_dt | delta_id >= 1);
  num_passes = length(I_break) + 1;
  pass_numbers = ones(size(gps_time,1),1);

  % Get the start and stop indices for each pass
  start = [1; I_break+1];
  stop = [I_break; size(gps_time,1)];

  if nargout > 1,
    % Get the start and stop times and id numbers
    pass_times = [[1:num_passes]' sec2gpst(data(start,1)) data(stop,1)-data(start,1) ...
                        data(start,3:4)];
    if nargout == 3 & nargin == 4,  % Compute the pass_summary from other_vis_data information
       % other_vis_data values at start and stop of passes
       pass_summary(:,:,1) = [data(start,5:4+size(other_vis_data,2))];
       pass_summary(:,:,2) = [data(stop,5:4+size(other_vis_data,2))];
       % Fill in min and max values for first pass
       pass_summary(1,:,3) = [min(data(start(1):stop(1),5:4+size(other_vis_data,2)))];
       pass_summary(1,:,4) = [max(data(start(1):stop(1),5:4+size(other_vis_data,2)))];
    end;
  end;
  for j=2:num_passes-1,
%    pass_numbers(I_break(j-1)+1:I_break(j)) = j;
     pass_numbers(start(j):stop(j)) = j;
     if nargout == 3 & nargin == 4,  % Compute the pass_summary from other_vis_data information
       % Fill in min and max values for each pass pass
       pass_summary(j,:,3) = [min(data(start(j):stop(j),5:4+size(other_vis_data,2)))];
       pass_summary(j,:,4) = [max(data(start(j):stop(j),5:4+size(other_vis_data,2)))];
     end;
  end;
  if num_passes > 1,
%    pass_numbers(I_break(num_passes-1)+1:size(gps_time,1)) = num_passes;
     pass_numbers(start(num_passes):stop(num_passes)) = num_passes;
     if nargout == 3 & nargin == 4,  % Compute the pass_summary from other_vis_data information
       % Fill in min and max values for each pass pass
       j = length(start);
       pass_summary(j,:,3) = [min(data(start(j):stop(j),5:4+size(other_vis_data,2)))];
       pass_summary(j,:,4) = [max(data(start(j):stop(j),5:4+size(other_vis_data,2)))];
     end;
  end;
else,  % Only one pass
  pass_numbers = ones(size(gps_time,1),1);
  if nargout > 1,
    pass_times(1,:) = [1 sec2gpst(min(data(:,1))) max(data(:,1))-min(data(:,1)) 1 1];
    if nargout == 3 & nargin == 4,  % Compute the pass_summary from other_vis_data information
       % other_vis_data values at start and stop of passes
       pass_summary(:,:,1) = [data(1,5:4+size(other_vis_data,2))];
       pass_summary(:,:,2) = [data(size(data,1),5:4+size(other_vis_data,2))];
       % Fill in min and max values for first pass
       pass_summary(:,:,3) = [min(data(:,5:4+size(other_vis_data,2)))];
       pass_summary(:,:,4) = [max(data(:,5:4+size(other_vis_data,2)))];
    end;
  end;
end;

% Reorder output data to match input
pass_numbers(i_2) = pass_numbers;
pass_numbers(i_1) = pass_numbers;

% Sort the pass_times and pass_summary data according to start times
if nargout > 1,
  [pass_sec, I_pass] = sort(gpst2sec(pass_times(:,2:3)));
  pass_times = pass_times(I_pass,:);
  if nargout==3 & nargin==4,
    pass_summary = pass_summary(I_pass,:,:);
  end;
end;

%%%% END ALGORITHM CODE %%%%%

% end PASSDATA
