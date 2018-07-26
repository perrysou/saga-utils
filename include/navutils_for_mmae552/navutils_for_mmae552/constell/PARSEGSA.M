function gsa_out = parsegsa(line)

% gsa_out = parsegsa(line);
%
% Function to parse GSA NMEA messages.
%
% Input:
%   line    - string with NMEA GSA data for one message
% Output:
%   gsa_out - output data message for GSA data type
%               [3, gsa_data] 3 is the internal GSA data type
%              where gsa_data is an nx17 matrix with columns
%               [auto, fix_dim, prn_1, prn_2, ..., prn_12, PDOP, HDOP, VDOP]
%              auto - 0 = auto 2/3-D mode, 1 = manual 2/3-D mode, -1 = unknown
%              fix_dim is the dimension fix for this data, 2 = 2-D, 3 = 3-D 
% 
% Note: prn values that are not filled in the GSA message will be filled 
%       with -1. PDOP values for 2D fixes will also be filled with -1. When
%       the fix is not 2D or 3D, the HDOP and VDOP will be filled with -1.
%                                      
% See also READNMEA, PARSNMEA, PARSEGGA, PARSEGSV 

% Written by: Jimmy LaMance 2/18/97
% Copyright (c) 1998 by Constell, Inc.

% functions called: NMEANEXT

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gsa_out=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PARSEGSA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% verify that line is a string
if ~isstr(line)
  fprintf('NMEA line input to PARSEGSA must be a string. \n');
  fprintf('See help PARSEGSA for details.') 
  if DEBUG_MODE
    fprintf('Error from PARSEGSA:  ')
    fprintf('Wrong type of line variable to PARSEGSA.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Wrong type of line variable to PARSEGSA.');
  end % if DEBUG_MODE
end % if ~isstr(line)

%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

ll = size(line,2);     % length of line
gsa_out = [];          % initialize the output to be blank

% default size of message output
gsa_size = 18;     % size of GSA out data is 18, checksum not supported

if strcmp(line(1:6),'$GPGSA') == 1    % found a GSA message 
  % allocate parse_out to be the right dimension
  gsa_out = ones(gsa_size,1) * inf;
  gsa_out(1) = 3;         % message type 1, GSA 

  [Auto, index] = nmeanext(line, 8); 
  if strcmp(Auto,'A')
    gsa_out(2) = 0;        % auto selection of 2D or 3D fix
  elseif strcmp(Auto,'M')
    gsa_out(2) = 1;        % manual selection of 2D or 3D fix
  else
    gsa_out(2) = -1;       % unknown selection of 2D or 3D fix
  end % if
  
  % fix dimension (3 -> 3-D, or 2 -> 2-D)
  [fix_dim, index] = nmeanext(line, index+1);    
  
  gsa_out(3) = fix_dim;
  
  for i = 4:15    % loop over the 12 places for satellites 
    [sv_num, index] = nmeanext(line,index+1);
    if isempty(sv_num) | sv_num == inf
      gsa_out(i) = -1;
    else
      gsa_out(i) = sv_num;
    end % if sv_num ~= []
  end % for i = 4:15    % loop over the 12 places for satellites 
    
  % PDOP 
  if gsa_out(3) ~= 3   % if it's not a 3-D fix, then there is no PDOP
    gsa_out(16) = -1;  % set the output data to the flag value
  else    
    [pdop, index] = nmeanext(line, index+1);   % read in the PDOP 
    if strcmp(pdop,'')                         % verify that there is a value
      gsa_out(16) = -1;
    else
      gsa_out(16) = pdop;                      % assign the value to the data
    end % if sv_num ~= []
  end % if gsa_out(3) ~= 3
  
  if gsa_out(3) ~= 3 &  gsa_out(3) ~= 2       % not a valid position fix
    gsa_out(17) = -1;       % flag the HDOP as invalid
    gsa_out(18) = -1;       % flag the VDOP as invalid
  else    
    % HDOP
    [gsa_out(17), index] = nmeanext(line, index+1);    
  
    % VDOP
    [gsa_out(18), index] = nmeanext(line, index+1);  
  end % if gsa_out(3) ~= 3 |  gsa_out(3) ~= 2      
  
else
  fprintf('Warning from PARSEGSA. GSA message not found in the current line.\n');
  fprintf('  %s\n\n',line);
  
end % if strcmp(line(1:6),'$GPGSA') == 1 

%%%%% END ALGORITHM CODE %%%%%

% end of PARSEGSA
  



  
