function gga_out = parsegga(line)

% gga_out = parsegga(line);
%
% Function to parse GGA NMEA messages.
% 
% Input:
%   line    - string with one NMEA GGA data message
% Output:
%   gga_out - output data message for GGA data type
%                        1 1  2  3   4   5   6      7        8      9      10  
%           columns are [1 hh mm ss lat lon alt diff_flag num_used HDOP dgps_age]
%             lat and lon are in deg, hgt and geoid_height are in m
%             diff_flag - 1 = GPS, 2 = DGPS 
% 

% Written by: Jimmy LaMance 4/7/97 
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gsa_out=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PARSEGGA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% verify that line is a string
if ~isstr(line)
  fprintf('NMEA line input to PARSEGGA must be a string. \n');
  fprintf('See help PARSEGGA for details.') 
  if DEBUG_MODE
    fprintf('Error from PARSEGGA:  ')
    fprintf('Wrong type of line variable to PARSEGGA.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Wrong type of line variable to PARSEGGA.');
  end % if DEBUG_MODE
end % if ~isstr(line)

%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

ll = size(line,2);     % length of line
gga_out = [];                 

% default size of message output
gga_size = 10;

if strcmp(line(1:6),'$GPGGA') == 1    % found an Ashtech GGA message 
  
  % Find all of the commas in the line
  I_comma = findstr(line,',');
  
  % A valid GGA message has 14 commas
  if length(I_comma) ~= 14
    return
  end
  
  % parse out the time
  time_field = str2num(line(I_comma(1)+1:I_comma(2)-1));
  gga_out(1) = fix(time_field / 10000);    % hours
  gga_out(2) = fix(time_field / 100) - gga_out(1) * 100;    % minutes
  gga_out(3) = time_field - gga_out(1) * 10000 - gga_out(2) * 100;   % seconds

  lat = str2num(line(I_comma(2)+1:I_comma(3)-1));
  North = line(I_comma(3)+1:I_comma(4)-1);
  deg = fix(lat / 100);
  min = lat - deg * 100;

  if strcmp(North,'N')
    gga_out(4) = deg + min / 60;      % latitude in degrees
  elseif strcmp(North,'S')
    gga_out(4) = -1 * (deg + min / 60);      % (south) latitude in degrees
  else
    gga_out(4) = inf;
  end % if
  
  lon = str2num(line(I_comma(4)+1:I_comma(5)-1)); 
  East = line(I_comma(5)+1:I_comma(6)-1);  
  deg = fix(lon / 100);
  min = lon - deg * 100;

  if strcmp(East,'E')
    gga_out(5) = deg + min / 60;      % longitude in degrees 
  elseif strcmp(East,'W')
    gga_out(5) = -1 * (deg + min / 60);      % (west) longitude in degrees 
  else
    gga_out(5) = inf;
  end % if

  % 0 - invalid, 1 - GPS, 2 - DGPS  
  
  valid = str2num(line(I_comma(6)+1:I_comma(7)-1));    % fix quality indicator 
  if ~isempty(valid)
    gga_out(7) = valid;
  else
    gga_out(7) = inf;
  end
  
  % number of satelltes
  num_sats = str2num(line(I_comma(7)+1:I_comma(8)-1));    % number of satellite 
  if ~isempty(num_sats)
    gga_out(8) = num_sats;
  else
    gga_out(8) = inf;
  end
  
  % HDOP
  hdop = str2num(line(I_comma(8)+1:I_comma(9)-1));    
  if ~isempty(hdop)
    gga_out(9) = hdop;
  else
    gga_out(9) = inf;
  end
  
  % height above the Ellipsoid
  hgt = str2num(line(I_comma(9)+1:I_comma(10)-1)); 
  if ~isempty(hgt)
    gga_out(6) = hgt;
  else
    gga_out(6) = inf;
  end
  
  % Geoid height, not saved to the output
%  dummy = line(I_comma(11)+1:I_comma(12)-1);
%  meter_flag = line(I_comma(12)+1:I_comma(13)-1); 

  % time since last DGPS update (blank if not using DGPS
  dgps_time = str2num(line(I_comma(13)+1:I_comma(14)-1));
  if isempty(dgps_time)
    gga_out(10) = inf;
  else
    gga_out(10) = dgps_time;
  end % if
  
  % DGPS Base station number (blank if not using DGPS), not decoded here, may 
  % be needed later when testing at Holloman with multiple base stations
  I_star = findstr(line,'*'); 
  
  if length(I_star) ~= 1   % too many start characters
    gga = [];
    return
  end
  
%  base_id = str2num(line(I_comma(14)+1:I_star-1));
%  if isempty(base_id)
%    gga_out(13) = inf;
%  else
%    gga_out(13) = base_id;
%  end % if
  
  % checksum
  check_sum = line(I_star+1:length(line));
  
  % get all of the characters between the header and the star before the check sum
  start_ind = 2;                    % exclude the $ sign
  stop_ind = length(line) - 3;      % remove the * and 2 digit hex check sum
  
  chk_line = line(start_ind:stop_ind);
  
  % Compute the check sum
  s_chk = double(chk_line(1));
  for i = 2:length(chk_line)
    s_chk = bitxor(double(chk_line(i)),s_chk);
  end      
  
  valid = 0;    % set this as an invalid message until the checksum has passed
  s_chk = dec2hex(s_chk);    % convert to hex
  
  if ~isstr(check_sum)
    check_sum = num2str(check_sum);
  end % if ~isstr(check_sum)

  if strcmp(check_sum,s_chk)
    valid = 1;
  else
    valid = 0;
  end % if strcmp(check_sum,s_chk)
    
  % make sure that at least valid time, lat and lon data were obtained, 
  % otherwise return a blank message
  if any(find(gga_out(1:9) == inf)) | valid == 0
%    fprintf('Error parsing GGA message\n  %s\n',line);
%    keyboard
    gga_out = [];
  end % if any(find(gga_out(1:9) == inf))
  
else
  fprintf('Warning from PARSEGGA. GGA message not found in the current line.\n')
  fprintf('  %s\n\n',line);

end % if strcmp(line(1:6),'$GPGGA') == 1   

% Add the GGA flag              
gga_out = [1 gga_out];

%%%%% END ALGORITHM CODE %%%%%

% end of PARSEGGA


  



  
