function stop_flag = err_chk(estruct);

% stop_flag = err_chk(estruct);
%
% Input variable error checking. Checks for matrix dimensions meeting
% required input dimension, common matrix dimensions between variables,
% and NaN, inf, and real checking.  All errors report messages to the screen.
% Dimension failures caues the stop_flag to be set.  NaN, inf, and real failues
% cause a warning message only.
% 
% Input: estruct - error checking structure
%   estruct.func_name = 'string';              function name
%   estruct.variable(i).name = 'string';       i-th variable name from input
%   estruct.variable(i).req_dim = nxm array;   i-th input required dimensions (opt.)
%   estruct.variable(i).var = jxk;             i-th variable
%   estruct.variable(i).type = 'type_string';  i-th variable type (optional)
%
%   The required dimensions (req_dim) is nxm where multiple options are allowed.
%   For example, a valid input could be 1x3 or 1x4 which would lead to the
%     req_dim = [1 3; 1 4];
%   For matching inputs between multiple variables, use a number greater than 900
%   to indicate that the dimensions must match.  For example if variable #1 and
%   variable #2 must match in the row dimension, but have different requirements
%   the column dimension, it would be handled as follows
%     estruct.variable(1).req_dim = [901 2];
%     estruct.variable(2).req_dim = [901 3; 901 4];
%   This input means that the variables #1 and #2 must match in the row dimension 
%   and the column dimensions are 2 and 3 or 4, respectively.  
%   There is no restriction on the number of variables, the number of variables
%   that must match a given dimension, or the number of allowed dimensions in 
%   the required dimension (req_dim) field.
%
%   The 'type' of variable input is used for the sanity checking of the values
%   of that variable.  Following are the current valid types and thier bounds.
%      LLA_RAD 
%          latitude bounds = [-pi/2 pi/2], longitude boundes = [-pi 2*pi] 
%      GPS_TIME
%          GPS week bounds = [0 3640], GPS sec bounds = [0 604800]
%      ANGLE_RAD
%          radian based angle bounds = [-2*pi 2*pi]
%      ELEVATION_RAD
%          elevation angle in radians bounds = [-pi/2 pi/2]
%      ELEVATION_RAD_0
%          elevation angle in radians bounds = [0 pi/2]
%      KEPLER_ORBIT
%          Keplerian orbit check, input [a e], (1x2) or (nx2)
%          a is in meters, e is dimensionless.  Checks perigee and apogee to be 
%          above the surface of the Earth.  If neither apsis is above the surface
%          of the Earth a warning message is issued.
%      STRING 
%          variable must be a string 
%
%   If an invalid or unknown type is specified a warning message will be issued.
%   The type checking is completely optional.
%
% Output: stop_flag - terminal condition flag, 1 = stop the called function
%                     otherwise continue execution  
% 
% Note: If the global variable DEBUG_MODE is set to 1, the stop_flag will
%       not be set and execution will continute.  Error messages will be printed
%       to the screen.

% Written by: Jimmy LaMance 4/15/98
% Copyright (c) 1998 by Constell, Inc.            

% functions called: none

global DEBUG_MODE           

% if the DEBUG_MODE has not been set, set it to 0, the default error checking mode
if isempty(DEBUG_MODE)
  DEBUG_MODE = 0;
end % if isempty(DEBUG_MODE)

% Find the number of variable that are being checked for correct dimensions
num_to_check = size(estruct.variable,2);

% Initialize the variable for similar dimension comparison between variables
sim_compare = [];
var_compare = [];
stop_flag = 0;

% Loop over the variables to check and flag bad inputs
for i = 1:num_to_check
  % Check to see if this variable has an input required dimension
  if isfield(estruct.variable(i),'req_dim')
    if ~isempty(estruct.variable(i).req_dim) 
      valid_dim = 0;     % flag this as invalid
      estruct.variable(i).valid = 0;
    
      % Get the input and required dimension in new variable names for ease of use
      estruct.variable(i).in_dim = size(estruct.variable(i).var);
      in_dim = estruct.variable(i).in_dim;
      req_dim = estruct.variable(i).req_dim;
    
      % Make sure that the input required dimension and input variable dimensions
      % match in array sizes (?x2)
      if size(req_dim,2) == size(in_dim,2)
    
        % Loop over the required dimension to see if there are any matches.  A valid
        % match at this point includes a match in a single dimension if the other 
        % dimensions are variable (900+).
        for jjj = 1:size(req_dim,1)                    
          % Check for an exact match (e.g. 1x3)
          if sum(abs(in_dim - req_dim(jjj,:))) == 0   
            valid_dim = 1;             
            estruct.variable(i).valid = 1;
          end

          % Check for a single dimension match with a 900 in the other dimension
          I_900 = find(req_dim(jjj,:) >= 900);
          I = find(req_dim(jjj,:) < 900);
          if ~isempty(I) & ~isempty(I_900)
            if sum(abs(in_dim(I) - req_dim(jjj,I))) == 0 
              valid_dim = 1;                     
            end % if sum(abs(in_dim(I) == req_dim(jjj,I))) == 0
          end % if in_dim(jjj,:) == req_dim(jjj,:) 

        end % for i = 1:size(req_dim,1)
      end % if size(req_dim,2) == size(in_dim,2)

      % If this is an invalid input, flag it as so
      if valid_dim == 0
       fprintf('%s input %s dimension is incorrect. \n',...
                 upper(estruct.func_name),estruct.variable(i).name);
        fprintf('Input was %d x %d.  Valid %s dimensions are\n',...
                 in_dim,estruct.variable(i).name);

        % Print valid dimensions to the screen
        for k = 1:size(estruct.variable(i).req_dim,1)
          valid_dim_num = estruct.variable(i).req_dim(k,:);
      
          % See if any of the valid dimensions have > 900 and substitute with n
          I_n = find(valid_dim_num > 900);
      
          % Replace the 900 with a single digit 9.  If the 900 stays then
          % the conversion to a string ends up with 3 digits.  This is the simple
          % way to avoid that.  The 9 is just a placeholder at this point.
          if ~isempty(I_n)
            valid_dim_num(I_n) = ones(size(I_n)) * 9;
          end
        
          % Convert the valid dimensions to a string
          valid_dim = num2str(valid_dim_num);
      
          % Fill in the variable length array dimensions with the letter 'n' ASCII 110
          if ~isempty(I_n)
            valid_dim(I_n) = char(110);
          end
      
          % Print this set of valid dimensions
          fprintf('    %s   =>   %c   x %s\n', ...
                   estruct.variable(i).name, deblank(valid_dim'));
        end % for k = 1:length(vars)

        fprintf('Type "help %s" for correct dimensions and details.\n\n', ...
                 lower(estruct.func_name));
    
        dbstack
        fprintf('\n');

        if DEBUG_MODE == 1
          fprintf('Error from %s:  ',upper(estruct.func_name))
          fprintf('Wrong size on %s variable to %s.\n',...
                       estruct.variable(i).name,upper(estruct.func_name));
        else
           % Set the flag to terminate the function upon return
          stop_flag = 1;
        end % if DEBUG_MODE
      end % if size(t2,2) ~= 2

      % Now mark this variable if it has to be size compared to another variable
      I_sim = find(estruct.variable(i).req_dim > 900);          
      if ~isempty(I_sim) & estruct.variable(i).valid ~= 1 & valid_dim == 1
        sim_compare = [sim_compare; estruct.variable(i).req_dim(I_sim)]; 
        var_compare = [var_compare; ones(size(I_sim)) * i];
      end                         
    end % if isfield(estruct.variable(i),'req_dim') 
  end % isfield(estruct.variable(i),'req_dim')
end % for i = 1:num_to_check

% Compare similar dimension variables as flagged above.  If there are none,
% skip this section of the code.
if ~isempty(sim_compare)

  % Find all of the unique similar variable comparisons
  u_compares = unique(sim_compare);

  num_unique_compare = length(u_compares);

  for i = 1:num_unique_compare
    % Find the variables that need to be compared
    I_var = find(sim_compare == u_compares(i));
  
    vars = var_compare(I_var);
  
    % Get the input and required dimension for this set of compares
    in_dim = [];
    req_dim = [];
    for j = 1:length(vars)
      in_d_temp = ones(size(estruct.variable(vars(j)).req_dim,1), ...
                       size(estruct.variable(vars(j)).in_dim,1)) * ...
                       estruct.variable(vars(j)).in_dim;
      in_dim = [in_dim; in_d_temp]; 
      req_dim = [req_dim; estruct.variable(vars(j)).req_dim];
    end % for j = 1:length(vars)  
  
    % Find all of the required simensions < 900 (these don't have to match)
    I_match = find(req_dim > 900);
    if ~isempty(I_match)
      diffs = diff(in_dim(I_match));
    end % if ~isempty(I_match)

    if any(diffs ~= 0)
      fprintf('Dimensions on the following variables do not agree.\n')
      fprintf('Type "help %s" for correct dimensions and details.\n', ...
                 lower(estruct.func_name));
      unique_vars = unique(vars);
      for k = 1:length(unique_vars)
        fprintf('    %s    %d x %d\n',estruct.variable(unique_vars(k)).name, ...
                                   estruct.variable(unique_vars(k)).in_dim);
      end % for k = 1:length(vars)
      fprintf('\n');

      dbstack
      fprintf('\n');

      % Set the flag for stopping the calling function upon return, if the 
      % DEBUG_MODE is the default of 0
      if DEBUG_MODE == 0
        stop_flag = 1;
      end % if ~DEBUG_MODE                                        
  
    end % if any(diffs ~= 0) 
  
  end % for i = 1:num_unique_compare  
end % if ~isempty(sim_compare)
  
empty_flag = 0;
for i = 1:num_to_check
  if isempty(estruct.variable(i).var)                                 
    if DEBUG_MODE == 0 | DEBUG_MODE == 1
      fprintf('Warning from Constellation Toolbox Error Checking.\n');
      fprintf(' Empty variables found in input %s to function %s.\n', ...
               estruct.variable(i).name, upper(estruct.func_name));
      fprintf('Resulting use of these variables may cause subsequent ');
      fprintf('functions to fail.\n\n');

      dbstack
      fprintf('\n');
    end % if DEBUG_MODE == 0 | DEBUG_MODE == 1
     
    % Set the flag indicating that a variable is empty
    empty_flag = 1;
  end % if ~isempty(I_nan)
end % for i = 1:num_to_check

% See if there are any errors up to this point.  If there are, return
% to the calling function.
if empty_flag == 1 | stop_flag == 1
  return
end

% Check the input variables for reasonableness, NaN, inf, and real values
for i = 1:num_to_check                 
  % Start sanity check of variables that have a type supplied
  if isfield(estruct.variable(i),'type')                      
    if strcmp(estruct.variable(i).type,'LLA_RAD');
      % Check that the inputs are within the error bounds
      lla_rad_min = [-pi/2 -pi];     % min lat and lon
      lla_rad_max = [pi/2 2*pi];     % max lat and lon
      
      % Find out of bounds latitude values
      I_lat = find(estruct.variable(i).var(:,1) < lla_rad_min(1) | ...
                   estruct.variable(i).var(:,1) > lla_rad_max(1));
      
      if ~isempty(I_lat)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d latitude variable(s) out of %d\n',length(I_lat),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between -pi/2 and pi/2.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          dbstack
          fprintf('\n');
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)

      % Find out of bounds longitude values
      I_long = find(estruct.variable(i).var(:,2) < lla_rad_min(2) | ...
                   estruct.variable(i).var(:,2) > lla_rad_max(2));
      
      if ~isempty(I_long)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d longitude variable(s) out of %d\n',length(I_long),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between -pi and 2*pi.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('functions to fail.\n\n');
          dbstack
          fprintf('\n');
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)
                   
    elseif strcmp(estruct.variable(i).type,'GPS_TIME')
      % Check that the inputs are within the error bounds
      gps_time_min = [0 0];             % min GPS week and seconds
      gps_time_max = [3640 604800];     % max GPS week and seconds
      
      % Find out of bounds latitude values
      I_week = find(estruct.variable(i).var(:,1) < gps_time_min(1) | ...
                    estruct.variable(i).var(:,1) > gps_time_max(1));
      
      if ~isempty(I_week)      
        if DEBUG_MODE == 0
          fprintf('Error from Constellation Toolbox Error Checking.\n');
          fprintf('   %d GPS week variable(s) out of %d\n',length(I_week),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between 0 and 3640.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          
          dbstack
          fprintf('\n');
          % Invalid GPS times set the stop flag when in DEBUG_MODE == 0
          stop_flag = 1;
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_week)

      % Find out of bounds longitude values
      I_sec = find(estruct.variable(i).var(:,2) < gps_time_min(2) | ...
                   estruct.variable(i).var(:,2) > gps_time_max(2));
      
      if ~isempty(I_sec)      
        if DEBUG_MODE == 0        
          fprintf('Error from Constellation Toolbox Error Checking.\n');
          fprintf('   %d GPS seconds variable(s) out of %d\n',length(I_sec),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between 0 and 604800.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');

          dbstack
          fprintf('\n');
          % Invalid GPS times set the stop flag when in DEBUG_MODE == 0
          stop_flag = 1;
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)
    
    elseif strcmp(estruct.variable(i).type,'ANGLE_RAD')
      % Check that the inputs are within the error bounds
      ang_rad_min = -2*pi;             % min mod 2*pi angle
      ang_rad_max = 2*pi;              % max mod 2*pi angle
      
      % Find out of bounds latitude values
      I_ang = find(estruct.variable(i).var < ang_rad_min(1) | ...
                    estruct.variable(i).var > ang_rad_max(1));
      
      if ~isempty(I_ang)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d angular variable(s) out of %d\n',length(I_ang),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between -2*pi and 2*pi.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          
          dbstack
          fprintf('\n');
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)

    elseif strcmp(estruct.variable(i).type,'ELEVATION_RAD')
      % Check that the inputs are within the error bounds
      ang_rad_min = -pi/2;             % min mod -pi/2 angle
      ang_rad_max = pi/2;              % max mod pi/2 angle
      
      % Find out of bounds latitude values
      I_ang = find(estruct.variable(i).var < ang_rad_min(1) | ...
                    estruct.variable(i).var > ang_rad_max(1));
      
      if ~isempty(I_ang)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d elevation variable(s) out of %d\n',length(I_ang),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between -pi/2 and pi/2.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          
          dbstack
          fprintf('\n');
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)
    
    elseif strcmp(estruct.variable(i).type,'ELEVATION_RAD_0')
      % Check that the inputs are within the error bounds
      ang_rad_min = 0;                 % min 0 angle
      ang_rad_max = pi/2;              % max mod pi/2 angle
      
      % Find out of bounds latitude values
      I_ang = find(estruct.variable(i).var < ang_rad_min(1) | ...
                    estruct.variable(i).var > ang_rad_max(1));
      
      if ~isempty(I_ang)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d elevation variable(s) out of %d\n',length(I_ang),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid values are between 0 and pi/2.\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          
          dbstack
          fprintf('\n');
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_lat)
    
    elseif strcmp(estruct.variable(i).type,'KEPLER_ORBIT')
      % Check that the inputs are within the error bounds
      earth_radius = 6378137.0;             % WGS-84 Earth radius meters
      
      % Compute apogee radius
      %                   a
      ra = estruct.variable(i).var(:,1) .* (1 + estruct.variable(i).var(:,2));
      rp = estruct.variable(i).var(:,1) .* (1 - estruct.variable(i).var(:,2));

      % Find out of bounds latitude values
      I_apsis = find(ra < earth_radius & rp < earth_radius);
      
      if ~isempty(I_apsis)      
        if DEBUG_MODE == 0
          fprintf('Warning from Constellation Toolbox Error Checking.\n');
          fprintf('   %d Keplerian orbits orbit(s) out of %d\n',length(I_apsis),...
                   size(estruct.variable(i).var,1));
          fprintf('   were found to be out of bounds in input %s to function %s.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Valid orbits have an apogee or perigee greater than\n');
          fprintf('   the radius of the Earth (6378137 m).  Check inputs for\n');
          fprintf('   correct units on semi-major axis (a).\n');
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');

          dbstack
          fprintf('\n');
          
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_apsis)
    
    elseif strcmp(estruct.variable(i).type,'STRING')
      
      if ~isstr(estruct.variable(i).var)      
        if DEBUG_MODE == 0
          fprintf('Error from Constellation Toolbox Error Checking.\n');
          fprintf('   Input %s to function %s must be a string.\n', ...
                   estruct.variable(i).name, upper(estruct.func_name));
          fprintf('   Type "help %s" for correct units and details.\n', ...
                    lower(estruct.func_name));
          fprintf('   Resulting use of these variables may cause subsequent ');
          fprintf('functions to fail.\n\n');
          
          dbstack
          fprintf('\n');
          % Invalid GPS times set the stop flag when in DEBUG_MODE == 0
          stop_flag = 1;
        end % if DEBUG_MODE == 0
      end % if ~isempty(I_week)
    
    elseif ~isempty(estruct.variable(i).type)
      fprintf('Unknowm variable type %s for input variable %s from function %s.\n',...
               estruct.variable(i).type, estruct.variable(i).name, ...
               upper(estruct.func_name));
      fprintf('Check the error checking inputs for this function.\n\n');
      dbstack
      fprintf('\n');
    end % if strcmp(estruct.variable(i).type,'LLA_RAD')  

    % If the variable type was declared, check for NaN, inf, and real
    % Check for NaN
    I_nan = find(isnan(estruct.variable(i).var) == 1);
    
    % Check for Nan, if there are any, issues a warning
    if ~isempty(I_nan)                                 
      if DEBUG_MODE == 0 | DEBUG_MODE == 1
        fprintf('Warning from Constellation Toolbox Error Checking.\n');
        fprintf(' %d NaN variables found in input %s to function %s.\n', ...
                 length(I_nan),estruct.variable(i).name, upper(estruct.func_name));
        fprintf('Resulting use of these variables may cause subsequent ');
        fprintf('functions to fail.\n\n');
        dbstack
        fprintf('\n');
      end % if DEBUG_MODE == 0 | DEBUG_MODE == 1
    end % if ~isempty(I_nan)

    % Check for inf
    I_inf = find(isinf(estruct.variable(i).var) == 1);
   
    % Check for inf, if there are any, issues a warning
    if ~isempty(I_inf)                                 
      if DEBUG_MODE == 0 | DEBUG_MODE == 1
        fprintf('Warning from Constellation Toolbox Error Checking.\n');
        fprintf(' %d inf variables found in input %s to function %s.\n', ...
                 length(I_inf),estruct.variable(i).name, upper(estruct.func_name));
        fprintf('Resulting use of these variables may cause subsequent ');
        fprintf('functions to fail.\n\n');
        dbstack
        fprintf('\n');
      end % if DEBUG_MODE == 0 | DEBUG_MODE == 1
    end % if ~isempty(I_nan)
     
    % Check for real, if there are any, issues a warning
    if ~isreal(estruct.variable(i).var)                                 
      if DEBUG_MODE == 0 | DEBUG_MODE == 1
        fprintf('Warning from Constellation Toolbox Error Checking.\n');
        fprintf('     %s variable in input to function %s is not real.\n', ...
                 estruct.variable(i).name, upper(estruct.func_name));
        fprintf('Resulting use of these variables may cause subsequent ');
        fprintf('functions to fail.\n\n');
        dbstack
        fprintf('\n');
      end % if DEBUG_MODE == 0 | DEBUG_MODE == 1
    end % if ~isempty(I_nan)
  end % if ~isempty(estruct.variable(i).type)                      
end % for i = 1:num_to_check

% end of ERR_CHK
