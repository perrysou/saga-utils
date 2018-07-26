function [data, index] = nmeanext(line, start_index);

% [data, index] = nmeanext(line, start_index);
%
% Function to get the next field in an NMEA message string.
%
% Inputs:
%   line        - a single line with an NMEA data (string)
%   start_index - index into the line to begin searching for the data
% Outputs:
%   data        - next number or character string contained in the NMEA message
%   index       - starting index to be used to retrieve the next piece of data
%
% See also EX_NMEA, READNMEA, PARSNMEA, PARSEGGA, PARSEGSA, PARSEGSV 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
data=[]; index=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on NMEANEXT for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'NMEANEXT';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'start_index';
  estruct.variable(1).req_dim = [1 1];
  estruct.variable(1).var = start_index;
  
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% declare a global variable to track error messaging so that multiple
% errors are not reported for the same invalid NMEA message line
global LAST_LINE

% stop flag   
end_of_data = 0;
count = 1;       
line_index = start_index; 
ll = size(line,2);     % line length    

if start_index > ll
  if ~strcmp(LAST_LINE,line)
    fprintf('Already reached end of line in NMEANEXT.\n');
    fprintf('Requested starting index was %d.  ');
    fprintf('The length of the line was %d.\n',start_index,ll);
    fprintf('Input line was\n   %s\n',line);   
    fprintf('This is probably an abnormal end of a message.\n\n');
    data = inf; 
    index = ll;            % set the index to the end of the line 
    LAST_LINE = line;
    return           
  else
    data = inf; 
    index = ll;            % set the index to the end of the line 
    return           
  end % if ~strcmp(LAST_LINE,line)
    
end % if start_index > ll

% put an aritifcial comma at the end of the line so the parsing algorithm
% can easily detect the last bit of data.  this is required to handle the
% inclusion of the checksum at the end of a line
line = [line ','];

while end_of_data == 0   

  % check to see if the first character is a comma or *
  if count == 1 & strcmp(line(line_index),',') | ...
     strcmp(line(line_index),'*')  
    start_index = start_index + 1;  % set the start index to the next character
  end % if    

  % see if we've reached a comma yet or checksum flag or end of the line
  % we've found the end of this number
  if count ~= 1 & (strcmp(line(line_index),',') | ...  
     strcmp(line(line_index),'*') | ...     % end of message (begin checksum)
     line_index > ll)                       % end of line
    end_of_data = 1;                                                        
    index = line_index - 1;       % set the output to the current line_index

    if isempty(str2num(line(start_index:line_index-1)))
      data = line(start_index:line_index-1);              % character data 
      
      % search for any $ in the middle of a message 
      I_dollar = findstr(data,'$');

      if ~isempty(I_dollar)
        % there was a $ in the string, this is bad data
        data = inf; 
        return
      end % if ~isempty(I_dollar)

      % verify that there are no blanks in the character data
      I_blank = findstr(line(start_index:line_index-1),' ');

      % if there are blanks replace them with zeros      
      if ~isempty(I_blank)
        data_string = deblank(line(start_index:line_index-1));
        
        
        % search for missing characters and replace them with 0
        I_blanks = find(data_string == ' ');
        
        if ~isempty(I_blanks)     
          data_string(I_blanks) = ones(size(I_blanks,1),1) * 48;
          data_string = setstr(data_string);
        end
        
        % check to see if it's a number now that the spaces have been removed
        if ~isempty(str2num(data_string))
          data = str2num(data_string);
        else
          data = inf;
        end % if ~isempty(str2num(data_string))
      end % if ~isempty(I_blank)

    else    % numeric data
      if line_index ~= start_index 
        % get rid of end of line blanks
        data_string = deblank(line(start_index:line_index-1));
        
        % search for missing characters and replace them with 0
        I_blanks = find(data_string == ' ');
        
        if ~isempty(I_blanks)     
          data_string(I_blanks) = ones(size(I_blanks,1),1) * '0';
        end
        
        data = str2num(data_string);   % numeric data 
        
        % verify that there is a single numeric value to be returned
        if length(data) > 1   
          data = inf;
        end % if length(data) < 1
      else
        data = inf;   % just 2 commas next to each other, fill with flag data
      end % if
    end
    
    % make sure that the data is not empty, otherwise the calling routine 
    % may not be able to accomodate it
    if isempty(data)
      data = inf;
    end % if isempty(data)
    
    return
  end % if
  
  line_index = line_index + 1;
  count = count + 1;
  
end % while  

%%%%% END ALGORITHM CODE %%%%%

% end of NMEANEXT
  
