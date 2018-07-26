function [gga_out, gsa_out, gsv_out] = readnmea(file_name)

% [gga_out, gsa_out, gsv_out] = readnmea(file_name);
%
% Function to read NMEA data from a file.
%
% Inputs:
%   file_name - name of file to read NMEA data from (1xn) (string)
% Outputs: 
%   gga_out - NMEA data from $GPGGA messages
%              gga_out is an nx10 matrix with columns
%                 [hh mm ss lat lon alt diff_flag num_used HDOP dgps_age]
%              lat and lon are in deg, hgt and geoid_height are in m
%              fix_quality - 1 = GPS, 2 = DGPS
%   gsa_out - NMEA data from $GPGSA messages
%              gsa_out is an nx17 matrix with columns
%              [auto, fix_dim, prn_1, prn_2, ..., prn_12, PDOP, HDOP, VDOP]
%              auto - 0 = auto 2/3-D mode, 1 = manual 2/3-D mode, -1 = unknown
%              fix_dim is the dimension fix for this data, 2 = 2-D, 3 = 3-D 
%   gsv_out - NMEA data from $GPGSV messages
%              gsv_out is an nx49 matrix with columns
%              [num_sats, prn_1, el_1, az_1, snr_1, prn_2, ..., snr_12]
%              num_sats - number of satellite being tracked
%              elevation (el) and azimuth (az) are in deg 
%              signal-to-noise ration in dB (definition varies with receiver)
%
% Note: READNMEA does not parse the last line of the file.
% Note on GGA data: Invalid or empty values will be filled with inf.
% Note on GSA data: prn values that are not filled in the GSA message will be  
%       filled with -1. PDOP values for 2D fixes will also be filled with -1. 
%       When the fix is not 2D or 3D, the HDOP and VDOP will be filled with -1.
% Note on GVA data: Any unavailable satellite data will be filled with inf.
%
% See also EX_NMEA, PARSNMEA, PARSEGGA, PARSEGSA, PARSEGSV 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: PARSNMEA

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gga_out=[]; gsa_out=[]; gsv_out=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on READNMEA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% check inputs 
if ~isstr(file_name)
  fprintf('NMEA file_name input to READNMEA must be a string. \n');
  fprintf('See help READNMEA for details.') 
  if DEBUG_MODE
    fprintf('Error from READNMEA:  ')
    fprintf('Wrong type of file_name variable to READNMEA.\n');
    % return to the calling function without filling in the output variables
    return
  else                                      
    error('Wrong type of file_name variable to READNMEA.');
  end % if DEBUG_MODE
end % if ~isstr(line)

% check the size of the string
if size(file_name,1 ) ~= 1
  fprintf('NMEA file_name input to READNMEA must be (1xn). \n');
  frpintf('Input was %d x %d.\n',size(file_name));
  fprintf('See help READNMEA for details.') 
  if DEBUG_MODE
    fprintf('Error from READNMEA:  ')
    fprintf('Wrong size of file_name variable to READNMEA.\n');
    % return to the calling function without filling in the output variables
    return
  else                                      
    error('Wrong size of file_name variable to READNMEA.');
  end % if DEBUG_MODE
end % if size(file_name,1 ) ~= 1
%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% clear the global variable that is set in neamnext for error message supression
if exist('LAST_LINE')
  clear global LAST_LINE
end % if exist('LAST_LINE')

% initialize the output data to be blank
gga_out = []; 
gsa_out = []; 
gsv_out = [];

% open the file
fid = fopen(file_name, 'r');

% verify that the file was opened correctly
if fid < 1
  fprintf('Unable to open file named %s in READNMEA.\n',file_name);
  fprintf('Verify the file name and path.\n')
  error('Unable to open the specified file in READNMEA.');
end % if fid < 1

end_of_file = 0;   % end of file flag (1 = end of file)
line_count = 1;    % counter for the number of lines read
num_nmea = 0;      % counter for number of nmea messages read
obs_num = 1;       % observation number
gga_num = 1;       % number of current gga message 
gsa_num = 1;       % number of current gsa message 
gsv_num = 1;       % number of current gsv message 

while end_of_file == 0         % read in the whole data file

  % get a line from the file
  line = fgetl(fid);

  % if the current line is empty, if this is the end of the file
  % it will be detected with the end of file pointed check that follows
  if isempty(line)
    line = ' ';
  end % if isempty(line)

  % verfiy that it's not an end of file and that the file pointer
  % is located at the end of the file.  processing the last line
  % of an NMEA data file generates a large number of warnings unless
  % it happens to be a complete message, usually it's not.
  if line == -1 | feof(fid) == 1      % end of file marker
    fprintf('Number of lines read = %d\n',line_count);
    fprintf('Number of NMEA messages processed = %d\n\n', num_nmea);
    return 
  end % if

  % call the nmea parsing routine to see if this is a recogonized message
  parse_out = parsnmea(line, fid);         
  
  % if the parse_out is not empty, then it's found a recogonized message
  if ~isempty(parse_out) 
    num_nmea = num_nmea + 1;
  
    % add the parse output to the appropriate NMEA data output
    if parse_out(1) == 1   % GGA message 
      % load the parsed data into the gga message data output
      % columns are [1 hh mm ss lat lon alt diff_flag num_used HDOP dgps_age]
      gga_out(gga_num,:) = parse_out(2:11); 
      
      gga_num = gga_num + 1;
      obs_num = obs_num + 1;

    elseif parse_out(1) == 3   % GSA message 
      % load the parsed data into the gsa message data output
      % format [auto_flag, 3-D_indicator, prns, (spaces for 12),  ...
      %        pdop, hdop, vdop] 
      gsa_out(gsa_num,:) = parse_out(2:18)'; 
      
      gsa_num = gsa_num + 1;
      obs_num = obs_num + 1;

    elseif parse_out(1) == 2   % GSV message 
      % load the parsed data into the gsa message data output
      %  format [hours, minutes, seconds, lat, lon, hgt, ...
      %        fix_quality, num_tracked, HDOP, geoid]
      gsv_out(gsv_num,:) = parse_out(2:length(parse_out))'; 
      
      gsv_num = gsv_num + 1;
      obs_num = obs_num + 1;

    end % if
  
  end % if    

  % increase the line counter
  line_count = line_count + 1;
  
end % while end_of_file == 0  

% clear the global variable that is set in neamnext for error message supression
if exist('LAST_LINE')
  clear global LAST_LINE
end % if exist('LAST_LINE')

%%%%% END ALGORITHM CODE %%%%%

% end of READNMEA
