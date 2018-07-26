function [leap_sec] = utc2leap(UTC_time)

% [leap_sec] = utc2leap(UTC_time);
%
% Determines the number of leap seconds of offset between GPS and UTC time.
%
% Input:  
%   UTC_time - matrix of the form [year month day hour minute second]
%               with 4-digit year (1980), nx6 matrix
% Output: 
%   leap_sec - leap seconds relating UTC to GPS time
%
% See also UTC2GPS, GPS2UTC

% Written by: Jimmy LaMance 10/9/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

% data files called: leapsecs.dat
%   leapsecs.dat has the form of [year month day hour minute second leapsec]
%   where the year month day hour minute second are UTC times, 
%   year is four digits

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
leap_sec=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on UTC2LEAP for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'UTC2LEAP';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'UTC_time';
  estruct.variable(1).req_dim = [901 6];
  estruct.variable(1).var = UTC_time;
 
  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE) 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
% verify that the leapsecs.dat file exist
if exist('leapsecs.dat') ~= 2
  fprintf('Unable to open the data file LEAPSECS.DAT.\n');
  fprintf('Verify that the file is in the Matlab path.\n');
  if DEBUG_MODE == 1
    fprintf('Error msg: Unable to open LEAPSECS.DAT.\n')
  else 
    error('Unable to open LEAPSECS.DAT.');
  end % if DEBUG_MODE == 1
end % if exist('leapsecs.dat') ~= 2

% load leap seconds data file (leapsecs.dat)
load leapsecs.dat

% verify that the leapsecs loaded has the correct number of columns
if size(leapsecs,2) ~= 7
  fprintf('Error in LEAPSECS.DAT file.\n');
  fprintf('There were %d columns in the file.  There should be 7 columns.\n',size(leapsecs,2));
  fprintf('Verify the Matlab path is correct and the correct version of LEAPSECS.DAT is in use.\n');
  if DEBUG_MODE == 1
    fprintf('Error msg: Invalid LEAPSECS.DAT file.\n')
  else 
    error('Invalid LEAPSECS.DAT file.');
  end % if DEBUG_MODE == 1
end % if size(leapsecs,2) ~= 7
    
n_leaps = size(leapsecs,1);   % number of leap seconds given in the data file

% search leapsecs for each of the time given in UTC_time. Only have to search based on 
% year and month since the leap seconds are added at the beginning of Jan and July.                     
for jjj = n_leaps:-1:1
  I_n = find(UTC_time(:,1) + UTC_time(:,2) / 12 < leapsecs(jjj,1) + leapsecs(jjj,2) / 12 );
  leap_sec(I_n) = ones(size(I_n,1),1) * leapsecs(jjj,7) - 1;
end % for jjj 

% find 
I_l = find(UTC_time(:,1) + UTC_time(:,2) / 12 >= leapsecs(n_leaps,1) + leapsecs(n_leaps,2) / 12 );
leap_sec(I_l) = ones(size(I_l,1),1) * leapsecs(n_leaps,7);

% Reshape leap_sec so that it is nx1
leap_sec = leap_sec';

%%%%% END ALGORITHM CODE %%%%%

% end of UTC2LEAP
