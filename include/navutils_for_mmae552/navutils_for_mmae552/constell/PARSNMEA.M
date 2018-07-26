function parse_out = parsnmea(line, fid)

% parse_out = parsnmea(line, fid);
%
% Function to parse a single line NMEA message.
%
% Input:
%   line   - string with NMEA data for one message
%   fid    - file handle to the data file (optional) (required if processing
%             a GSV message because this message can span up to 3 lines for
%             a complete message)
% Output:
%   parse_out - output data message, varies with message type
%                [message_type, message_data]
%                supported message types include GPGGA, GPGSV, GPGSA
%                message_type -> 1 = GPGGA , 2 = GPGSV, 3 = GPGSA
%                message data for each NMEA message supported is described below
%                                      
%   gga -> message_data is an nx12 matrix with columns
%              [hours, mins, secs, lat, lon, hgt, fix_qual, num_svs, HDOP, ...
%               geoid_height, DGPS_latency, DGPS_station_number]
%              lat and lon are in deg, hgt and geoid_height are in m
%              fix_quality - 1 = GPS, 2 = DGPS
%   gsa -> message_data is an nx17 matrix with columns
%              [auto, fix_dim, prn_1, prn_2, ..., prn_12, PDOP, HDOP, VDOP]
%              auto - 0 = auto 2/3-D mode, 1 = manual 2/3-D mode
%              fix_dim is the dimension fix for this data, 2 = 2-D, 3 = 3-D 
%   gsv -> message_data is an nx49 matrix with columns
%              [num_sats, prn_1, el_1, az_1, snr_1, prn_2, ..., snr_12]
%              num_sats - number of satellite being tracked
%              elevation (el) and azimuth (az) are in deg 
%              signal-to-noise ration in dB (definition varies with receiver)
%
% See also READNMEA, PARSEGGA, PARSEGSA, PARSEGSV 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: PARSEGGA, PARSEGSV, PARSEGSA

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PARSNMEA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% verify that line is a string
if ~isstr(line)
  fprintf('NMEA line input to PARSNMEA must be a string. \n');
  fprintf('See help PARSNMEA for details.') 
  if DEBUG_MODE
    fprintf('Error from PARSNMEA:  ')
    fprintf('Wrong type of line variable to PARSNMEA.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Wrong type of line variable to PARSNMEA.');
  end % if DEBUG_MODE
end % if ~isstr(line)

% verify that the file ID is valid, if provided
if nargin > 1
  testfile = fopen(fid);
    
  if isempty(testfile)
    fprintf('Invalid file ID (fid) given to PASRNMEA.\n');
    fprintf('No GSV data will be processed.\n');
    parse_out = [];               % fill output with a blank matrix
    return
  end % if isempty(testfile)
end % if nargin > 1

%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% check the length of the input line
if length(line) < 6             % there's not enough data to be a full message 
  parse_out = [];        % fill output with a blank matrix
  return         
end % if length(line) < 6
  
% parse messages
if strcmp(line(1:6),'$GPGGA')          % found a GGA message 
  parse_out = parsegga(line);  

elseif strcmp(line(1:6),'$GPGSV')      % found a GSV message 
  if nargin < 2
    parse_out = parsegsv(line);        % call parsegsv without the file ID
  else  
    parse_out = parsegsv(line, fid);
  end % if nargin < 2
   
elseif strcmp(line(1:6),'$GPGSA')      % found a GSA message 
  parse_out = parsegsa(line);  

else                                   % no recgonized message found
  parse_out = [];               % fill output with a blank matrix
  
end % if  

%%%%% END ALGORITHM CODE %%%%%

% end of PARSNMEA
  



  
