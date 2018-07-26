function [err_string] = writyuma(ephemeris,filename,over_write)

% [err_string] = writyuma(ephemeris,filename,over_write);
%
% Write YUMA formatted almanacs for any constellation.
%
% Inputs:
%   ephemeris  - ephemeris matrix for all satellites (nx24), with columns of
%                 [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                 i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                 GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                 Ephemeris parameters are from ICD-GPS-200 with the 
%                 exception of perigee_rate. 
%   filename   - name of alamanac file to write. (optional) 
%                 The file naming convention of gps###.alm and 
%                 glo###.alm is preferrable for GPS and GLONASS constellations. 
%                 User defined constellations may want to use usr###.alm.  
%                 Default = 'default.alm' If almanacs are stored in a diectory, 
%                 use the full (or partial) path names to write almanac files 
%                 into that directory.  For example, using a PC ...
%                      writyuma('c:\mywork\thisproj\almanac\user989.alm')
%                 will write the specificed almanac (user989.alm) to the
%                 specified directory (c:\mywork\thisproj\almanac).
%   over_write - flag indicating that the file should be overwritten without 
%                a warning message if it already exists (optional).  
%                over_write = 1 -> supress the warning message if the file exists
%                over_write = 0 -> display the warning message if the file exists
% Outputs:
%   err_string - error message issued if the write failed.  This string is empty
%                 if writing the YUMA almanac was successful
%
% See also ALM2GEPH, KEP2GEPH

% Note: Yuma almanac format is as follows...       
%       [sv_num, health, ecc, GPS_sec, inc, asc_node_rate, sqrt_a, asc_node,...
%        perigee, mean_anomaly, Af0, Af1, GPS_week] 
%       with units of rad, s, and m^.5

% Written by: Jimmy LaMance May 1998
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output
err_string = [];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on WRITYUMA for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Check for optional variables and fill in the default values if not provided
if nargin < 2
  filename = 'default.alm';
end % if nargin < 2

if nargin < 3
  over_write = 0;
end % if nargin < 3

estruct.func_name = 'WRITYUMA';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'ephemeris';
estruct.variable(1).req_dim = [901 24];
estruct.variable(1).var = ephemeris;

estruct.variable(2).name = 'filename';
estruct.variable(2).type = 'STRING';
estruct.variable(2).var = filename;

estruct.variable(3).name = 'over_write';
estruct.variable(3).req_dim = [1 1];
estruct.variable(3).var = over_write;

% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
err_string = [];

% Check to see if the file to be written already exists
p = which(filename);

% If the variable p is not empty, then the file was found in the Matlab
% path.  If this file is in the current working directory, then issues
% a dialog box asking if it's OK to overwrite the file.
if ~isempty(p)
  % Check if the file was found in the current path. Strip off the file name
  % and the trailing directory delimiter.
  pp = fileparts(p);    % gets just the path
  
  % Compare the path of this file to the path where the almanac file
  % will be written.  Start be getting the path from the input filename.
  input_path = fileparts(filename);
  
  % If the input filename does not contain a path, the file will be written
  % into the current directory
  if isempty(input_path)
    input_path = pwd;
  end % if isempty(input_path)
  
  if strcmp(input_path,pp) & over_write == 0
    % There is a file that already exists in this location.  Prompt the user
    % to overwrite or canel.
    dlg_dtring = sprintf('The file %s alread exists.  Do you want to overwrite it?',...
                          filename);
    button = questdlg(dlg_dtring,'Almanac file write warning','Yes','No','No');
    
    if strcmp(button,'No')
      err_string = sprintf('No almanac file %s written in WRITYUMA.  File already exists.',...
                            filename); 
      return
    end % if strcmp(button,'No')
  end % if strcmp(input_path,pp)

end % if ~isempty(p)

% Open almanac file for writing
fid = fopen(filename,'w');

if fid == 11
  err_string = sprintf('Unable to open file %s.',filename);
  return
end % if fid_all < 1

% Loop over the satellite in the ephemeris and write the data into
% the almanac file.
num_sats = size(ephemeris,1);

for i = 1:num_sats
  % Write the header for this satellite
  write_string = sprintf('**** Week %d almanac for PRN-%d ***********',...
                          ephemeris(i,19),ephemeris(i,1));
  fprintf(fid,'%s\n',write_string);                          

  % Write the ID line for this satellite
  write_string = sprintf('ID:                          %03d',...
                          ephemeris(i,1));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Health line for this satellite.  Assumed healthy.
  write_string = sprintf('Health:                      000');
  fprintf(fid,'%s\n',write_string);                          

  % Write the Eccentricity line for this satellite
  write_string = sprintf('Eccentricity:                %11.10E',...
                          ephemeris(i,4));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Time of Applicability line for this satellite
  write_string = sprintf('Time of Applicability(s):    %11.6f',...
                          ephemeris(i,17));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Orbital Inclination(rad) line for this satellite
  write_string = sprintf('Orbital Inclination(rad):    %12.10f',...
                          ephemeris(i,7));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Rate of Right Ascen(r/s) line for this satellite
  write_string = sprintf('Rate of Right Ascen(r/s):    %11.10E',...
                          ephemeris(i,9));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the SQRT(A) (m^1/2) line for this satellite
  write_string = sprintf('SQRT(A) (m^1/2):             %13.7f',...
                          ephemeris(i,5));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Right Ascen at TOA(rad) line for this satellite
  write_string = sprintf('Right Ascen at TOA(rad):     %11.10E',...
                          ephemeris(i,6));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Argument of Perigee(rad) line for this satellite
  write_string = sprintf('Argument of Perigee(rad):    %11.10E',...
                          ephemeris(i,8));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Mean Anom(rad) line for this satellite
  write_string = sprintf('Mean Anom(rad):              %11.10e',...
                          ephemeris(i,2));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Af0(s) line for this satellite
  write_string = sprintf('Af0(s):                      %11.10e',...
                          ephemeris(i,21));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Af1(s/s) line for this satellite
  write_string = sprintf('Af1(s/2):                    %11.10e',...
                          ephemeris(i,22));
  fprintf(fid,'%s\n',write_string);                          
 
  % Write the Week line for this satellite
  write_string = sprintf('week:                        %d',...
                          ephemeris(i,19));
  fprintf(fid,'%s\n\n',write_string);                          
 
end % for i = 1:num_sats

fclose(fid);  % close the file
  
%%%%% END ALGORITHM CODE %%%%%

% end of WRITYUMA
