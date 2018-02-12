function [gps_alm, glo_alm] = readyuma(filename,rollover_flag)

% [gps_alm, glo_alm] = readyuma(filename,rollover_flag);
%
% Read in YUMA formatted GPS, GLONASS, or user created almanacs.
%
% Inputs:
%   filename - name of almanac file to read or GPS week number. (optional) 
%               The files must have the naming convention of gps###.alm and 
%               glo###.alm (or gls####.alm) if using the week number to specify a 
%               file for GPSor GLONASS. If a filename instead of a GPS week number  
%               is specified, only the gps_alm matrix will be returned and it 
%               will contain the almanac data for the given file. If not provided
%               the FIND_ALM function will locate the most recent GPS and GLONASS
%               almanacs. This function can also be used to read user created
%               almanacs, by providing the full almanac name and extension. The data
%               is returned in the gps_alm field.  Almanac names are stored as
%               mod(1024) weeks to conform to YUMA standards.
%   rollover_flag - flag indicating the number of GPS week rollovers that have
%               occured (1-st rollover on August 22, 1999) (optional, default = 1).
%               GPS time is kept as total number of GPS weeks (weeks go beyond 1024)
%               The GPS week in the almanac will be replaced with the GPS_week + 1024
%               for a rollover flag of 1.
%               
% Outputs:
%   gps_alm  - GPS almanac data matrix for all the GPS satellites found 
%               in the specified almanac file. (nx13)
%   glo_alm  - GLONASS almanac data matrix for all the GLONASS satellites 
%               found in the specified almanac file.  Supported when an 
%               almanac week number is provided. (mx13) 
%
% Note: Yuma almanac format is as follows...       
%       [sv_num, health, ecc, GPS_sec, inc, asc_node_rate, sqrt_a, ...
%        long of asc. node at weekly epoch,...
%        perigee, mean_anomaly, Af0, Af1, GPS_week] 
%       with units of rad, s, and m^.5
%
% See also FIND_ALM, KEP2GEPH

% Written by: Jimmy LaMance 11/1/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: FIND_ALM

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

if nargin < 2
  rollover_flag = 1;
end % if

if nargin == 0 
  if nargout == 2
    [inputfile(1,:), inputfile(2,:)] = find_alm; 
  else
    [inputfile(1,:)] = find_alm; 
  end % if nargout == 2    
end % if nargin == 0
  
% Check inputs to determine if a filename or week number was provided
if nargin >= 1
  if isstr(filename)        % if the input is a string
    inputfile = filename;   % set the read file name to the input
  else                      % else, ...generate the file name based on input week
   
    % do a sanity check on the week number
    GPS_week_max = 3640;    % maximum value of GPS weeks
    GPS_week_min = 0;       % minimum value of GPS weeks
  
    if filename > GPS_week_max | filename < GPS_week_min
      fprintf('The GPS week provided (%d) is out of bounds. \n',filename);
      fprintf('Mimimum and maximum values are %d and %d.\n',...
               GPS_week_max,GPS_week_min);
      if DEBUG_MODE
        fprintf('Error message from READYUMA: \n');
        fprintf('The GPS week variable (filename) to READYUMA is out of bounds.\n');
        return 
      else
        error('The GPS week variable (filename) to READYUMA is out of bounds.')
      end % if DEBUG_MODE
    end % if filename > GPS_week_max | filename < GPS_week_min 

  
    % set up the file name for the GPS almanac using the user provided week
    inputfile(1,:) = sprintf('gps%d.alm',filename);   % GPS almanac
 
    % if 2 output variables were requested, set up the file name for the 
    % GLONASS almanac also
    if nargout == 2
      inputfile(2,:) = sprintf('glo%d.alm',filename);   % GLONASS almanac
    
    end % if nargout == 2
  end % if isstr(filename)
end % if nargin == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% compute how many almanac files there are based on the size of the inputfile
num_almanac_files = size(inputfile,1);

% Open almanac file(s) for reading
for i = 1:num_almanac_files
  fid_all(i) = fopen(inputfile(i,:),'r');
end % for i = 1:num_almanac_files  

% Verify that the opens were successful and set up which files will be read
if num_almanac_files == 1
  
  % check for a valid file ID for the GPS almanac
  if fid_all(1) < 0     
    fprintf('Unable to open the almanac file %s. \n', inputfile(1,:))
    fprintf('Check that the file name and path are correct.\n');
    fprintf('No GPS almanac will be returned\n')
    gps_alm = [];
    if DEBUG_MODE
      fprintf('Error message from READYUMA: \n');
      fprintf('Incorrect or invalid file/path name for GPS almanac.\n');
      return 
    else
      error('Incorrect or invalid file/path name for GPS almanac.')
    end % if DEBUG_MODE
  end % if fid_all(1) < 0

% else if there are GPS and GLONASS almanac files
elseif num_almanac_files == 2
  
  % check that both GPS and GLONASS almanacs were opened correctly
  if fid_all(1) < 0 & fid_all(2) < 0
    fprintf('Unable to open a GPS or GLONASS almanac file %s %s. \n', inputfile(1,:),inputfile(2,:))
    fprintf('Check that the file names and path are correct.\n');
    fprintf('No almanac data will be returned\n\n')
    gps_alm = [];
    glo_alm = [];
    if DEBUG_MODE
      fprintf('Error message from READYUMA: \n');
      fprintf('Incorrect or invalid file/path name for GPS and GLONASS almanacs.\n');
      return 
    else
      error('Incorrect or invalid file/path name for GPS and GLONASS almanacs.')
    end % if DEBUG_MODE
  
  % check that the GPS almanac was opened correctly
  elseif fid_all(1) < 0
    fprintf('Unable to open the almanac file %s. \n', inputfile(1,:))
    fprintf('Check that the file name and path are correct.\n');
    fprintf('No GPS almanac will be returned.\n\n')
    gps_alm = [];

  % check that the GLONASS almanac was opened correctly
  elseif fid_all(2) < 0
    fprintf('Unable to open the almanac file %s. \n', inputfile(2,:))
    fprintf('Check that the file name and path are correct.\n');
    fprintf('No GLONASS almanac will be returned.\n\n')
    glo_alm = [];

  end % if fid_all(1) < 0 & fid_all(2) < 0

end % if num_almanac_file == 1  

% set up a matrix of almanac files to read (don't read the ones with invalid
% file IDs)
I_valid_file = find(fid_all ~= -1);
alm_files_to_read = length(I_valid_file);

for k = 1:alm_files_to_read
  fid = fid_all(I_valid_file(k));     % set the file identifier for this loop
  
  % Read through the file and fill in the almanac matrix

  err_msg = blanks(0);                 % set the error message to blank  
  end_of_file = 0;              % set the end of file flag to 0
                                % using the not makes the logic a little easier 
  read_error = 0;               % set a flag if a read error has occurred 
  num_lines_per_almanac = 13;   % 15 lines for each sat. in the almanac file
                              
  satellite_count = 1;          % counter with the the current satellite count
  almanac_head = fgetl(fid);    % get the header line for the first satellite

  if ~isstr(almanac_head)       % does the line read contain a string
    if almanac_head == -1       % is there an end of file on the first line
      end_of_file = 1;
      read_error = 1;           % read error in the given almanac file
    end % if almanac_head == -1
  end % if ~isstr(almanac_head)
  
  clear alm                     % clear any previous almanacs
  
  while end_of_file == 0 & read_error == 0
  
    for i = 1:num_lines_per_almanac
      almanac_line = fgetl(fid);      % get the next line from the almanac file

      if almanac_line == -1           % is there an end of file
        end_of_file = 1;              % set end of file flag

      else
        I_colon = find(almanac_line == ':');  % find the colon in the line
        if isempty(I_colon)   				  % search for a '('
          I_colon = find(almanac_line == ')');
        end
        
        I_end = size(almanac_line,2);         % find the end of the line
      
        if any(I_colon) & I_end > 3              % valid almanac line
           [alm(satellite_count,i),read_count,err_msg] = ...
            sscanf(almanac_line(I_colon+1:I_end),'%f');    
        else
          read_error = 1;                     % set the read error flag
          keyboard
        end % if any(I_colon) & I_end > 3 

        if ~isempty(err_msg)                 % verify the read was successfull
          read_error = 1;                   % set no_read_error flag
          i = num_lines_per_almanace + 1;   % set loop parameter to exit loop
        end % if any(err_msg)

      end % if almanac_line == -1                                            
    
    end % for i = 1:num_lines_per_almanac                                                                        

    almanac_tail = fgetl(fid);       % get the tailer line for this satellite
    almanac_head = fgetl(fid);       % get the header line for next satellite

    if ~isstr(almanac_tail) | ~isstr(almanac_head)
      if almanac_head == -1 | almanac_tail == -1      % is there an end of file 
        end_of_file = 1;
      end % if almanac_head == -1
    end % if ~isstr(almanac_tail) | ~isstr(almanac_head) 
  
    
    satellite_count = satellite_count + 1;  % increment the satellite counter
    
  end % while end_of_file == 0 & no_read_error  

  if I_valid_file(k) == 1      % the almanac should fo in the GPS almanac return
    gps_alm = alm;
  
  elseif I_valid_file(k) == 2  % the almanac should fo in the GLONASS almanac return
    glo_alm = alm;
  
  end % if k == 1
  
  fclose(fid);  % close the file just read
  
end % for k = 1:num_almanac_files

if rollover_flag ~= 0
  gps_alm(:,13) = gps_alm(:,13) + rollover_flag*1024;
  glo_alm(:,13) = gps_alm(:,13) + rollover_flag*1024;
end % if rollover_flag ~= 0

%%%%% END ALGORITHM CODE %%%%%

% end of READYUMA
