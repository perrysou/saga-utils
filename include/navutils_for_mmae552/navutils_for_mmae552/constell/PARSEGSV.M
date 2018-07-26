function gsv_out = parsegsv(line, fid)

% gsv_out = parsegsv(line, fid);
%
% Function to parse GSV NMEA messages.
%  
% Input:
%   line    - string with NMEA GSV data for one message
%   fid     - file handle to the data file (optional) (required for processing
%              GSV messages that span more than a single line) 
% Output:
%   gsv_out - output data message for GSV message type
%              [2, gsv_data]   2 is the internal GSV message type
%              where gsv_data is an nx49 matrix with columns
%              [num_sats, prn_1, el_1, az_1, snr_1, prn_2, ..., snr_12]
%              num_sats - number of satellite in view
%              elevation (el) and azimuth (az) are in deg 
%              signal-to-noise ration in dB (definition varies with receiver)
% 
% Note: Any unavailable satellite data will be filled with inf.
%
% See also READNMEA, PARSNMEA, PARSEGSA, PARSEGGA 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: NMEANEXT

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gsv_out=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PARSEGSV for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% verify that line is a string
if ~isstr(line)
  fprintf('NMEA line input to PARSEGSV must be a string. \n');
  fprintf('See help PARSEGSV for details.') 
  if DEBUG_MODE
    fprintf('Error from PARSEGSV:  ')
    fprintf('Wrong type of line variable to PARSEGSV.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Wrong type of line variable to PARSEGSV.');
  end % if DEBUG_MODE
end % if ~isstr(line)

% verify that the file ID is valid, if provided
if nargin > 1
  testfile = fopen(fid);
    
  if isempty(testfile)
    fprintf('Invalid file ID (fid) given to PASREGSV.\n');
    fprintf('No GSV data will be processed.\n');
    gsv_out = [];
    return
  end % if isempty(testfile)

else
  fid = -1;             % set the file to an invalid file ID
end % if nargin > 1

%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% note on GSV messages...
% The GSV message spans multiple lines.  This routine is designed to have the
% first line (message) passed in with the file ID (obtained with fopen).  
% After determining that the line passed in is the first line of a GSV message,
% PARSEGSV will read in consecutive lines searching for the remaining set of
% GSV messages (there may be 1, 2, or 3 messages).  This is a variable
% length and variable number of lines message. 

% default size of message output
gsv_size = 50;     % size of the output gsv message

% save a copy of the first GSV message line.  
% used to for warnings and error messages
first_line = line;

if strcmp(line(1:6),'$GPGSV') == 1    % found a GSV message 
  % allocate parse_out to be the right dimension
  gsv_out = ones(gsv_size,1) * inf;
  gsv_out(1) = 2;         % message type 2, GSV 

  % number of lines for a full message
  [num_lines, index] = nmeanext(line, 8); 

  % the first number of lines, used to validity checking
  first_line_count = num_lines;
  
  % which line of the full message is this
  [sentence_number, index] = nmeanext(line, index + 1);
  
  % verify that this is the first line, otherwise return nothing
  if sentence_number ~= 1
    fprintf('Incorrect GPGSV sentence number found in PARSEGSV.\n');
    fprintf('  %s\n',line);                     
    fprintf('The first GSV sentence must be labeled as sentence number 1.\n');
    fprintf('This message is labeled message number %d.\n',sentence_number);
    fprintf('Check the input NMEA file (or string) for missing ');
    fprintf('GSV messages.\n\n');
    gsv_out = []; 
    return
  end % if line_number ~= 1  
  
  % the last sentence read, used for validity checking
  last_sentence = sentence_number;    
  
  % number of satellites to expect in the message
  [num_svs, index] = nmeanext(line, index + 1);
  gsv_out(2) = num_svs;
  
  % the first GSV message number of svs, used for validity checking 
  % this will be compared with the number of svs listed in later line
  % in the same GSV message
  first_sv_count = num_svs;

  % fill in the satellite data and read new lines as required
  for sv_count = 1:num_svs
    if sv_count == 5 | sv_count == 9   % need a new line from the file
      clear line    % clear the last line
      
      % verify that the file ID is valid (not set to -1), see error checking
      if fid == -1
        fprintf('No valid file was provided to read the remainder of the GSV ');
        fprintf('message.\n')
        fprintf('  %s\n',first_line);
        fprintf('There should be %d lines comprising this GSV message.\n',...
                 num_lines);
        fprintf('No GSV data will be returned for this message.\n\n');
        gsv_out = [];
        return
      end % if fid == -1
      
      line = fgetl(fid);
      
      % verify that the current line is not empty
      if isempty(line)
        fprintf('No line was readable from file %s.\n',fopen(fid));
        fprintf('Can not finish reading this GSV message.\n');
        fprintf('  %s\n',first_line);
        fprintf('Verify file name, path, and $GPGSV message formats.\n\n');
        gsv_out = []; 
        return 
      end % if isempty(line)
       
      % verfiy that it's not and end of file
      if line == -1 | feof(fid) == 1      % end of file marker
        end_of_file = 1;
        fprintf('\nEnd of file encountered in PARSEGSV.\n');
        fprintf('The last GSV message will not be processed.\n\n');
        gsv_out = []; 
        return 
      end % if
      
      % verfiy that the line is at least 6 characters long
      if length(line) < 6      % invalid
        fprintf('\nNMEA message shorter than 6 characters found in PARSEGSV.\n');
        fprintf('%s \n\n',line);
        fprintf('Unable to complete the read of this GSV message.\n');
        fprintf('  %s\n\n',first_line);
        gsv_out = []; 
        return 
      end % if
      
      % verify that this line is also a GSV message
      if ~strcmp(line(1:6),'$GPGSV') 
        % the Garmin 12 XL puts out PSLIB messages between the GSV messages
        % this check looks for the PSLIB message and skips over it before
        % determing that the GSV message is bad.  The PSLIB is a control
        % sentence that is used to tune beacon receivers
        while strcmp(line(1:6),'$PSLIB')
          clear line    % clear the last line
          line = fgetl(fid);
        end % while strcmp(line(1:6),'$PSLIB')
        
        % now that the PSLIB message have been passed, check again for a GSV
        if ~strcmp(line(1:6),'$GPGSV') 
          fprintf('Abnormal end to a GPGSV message. ');
          fprintf('No GSV header found.\n')
          fprintf('Unable to complete the read of this GSV message.\n');
          fprintf('  %s\n\n',first_line);
          gsv_out = [];   % return nothing
          return
        end % if ~strcmp(line(1:6),'$GPGSV')
        
      end % if strcmp(line(1:6),'$GPGSV') ~= 1 
      
      % begin parsing this line
      % number of lines for a full message
      [num_lines, index] = nmeanext(line, 8); 
      
      if first_line_count ~= num_lines
        fprintf('Abnormal end to a GPGSV message. ');
        fprintf('%s',line); 
        fprintf('The first message line was\n');
        fprintf('  %s\n',first_line);        
        fprintf('Incorrect total line number.\n\n')
        gsv_out = [];   % return nothing 
        return
      end % if first_line_count ~= num_lines
  
      % which line of the full message is this
      [sentence_number, index] = nmeanext(line, index + 1);
      
      if (sentence_number ~= last_sentence + 1)
        fprintf('Abnormal end to a GPGSV message. ');
        fprintf('%s',line);
        fprintf('Incorrect sentence number.  This should be sentence %d.\n\n',...
                 last_sentence + 1)
        gsv_out = [];   % return nothing
        return
      end % if (sentence_number ~= last_sentence + 1)      

      last_sentence = sentence_number;    % update the last_sentence variable

      % number of satellite to expect in the message
      [num_svs, index] = nmeanext(line, index + 1);
      
      if num_svs ~= first_sv_count
        fprintf('Abnormal end to a GPGSV message. ');
        fprintf('%s',line);
        fprintf('Incorrect satellite count.  ');
        fprintf('The satellite count should be %d.\n\n',num_svs)
        gsv_out = [];   % return nothing 
        return
      end % if num_svs ~= first_sv_count

    end % if sv_count == 5 | sv_count == 9

    % read in data for sv positions, sat number and snr
    gsv_index = (sv_count - 1) * 4 + 3;
    [gsv_out(gsv_index + 0), index] = nmeanext(line, index+1);     % prn
    [gsv_out(gsv_index + 1), index] = nmeanext(line, index+1);     % el (deg)
    [gsv_out(gsv_index + 2), index] = nmeanext(line, index+1);     % az (deg)
    [gsv_out(gsv_index + 3), index] = nmeanext(line, index+1);     % SNR

  end % for sv_count = 1:num_svs    
  
else
  fprintf('Warning from PARSEGSV. GSV message not found in the current line.\n')
  fprintf('  %s\n\n',line);
  
end % if strcmp(line(1:6),'$GPGSV') == 1 

%%%%% END ALGORITHM CODE %%%%%

% end of PARSEGSV

  



  
