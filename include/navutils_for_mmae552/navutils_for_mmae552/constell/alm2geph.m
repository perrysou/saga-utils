function [gps_ephem] = alm2geph(alm)

% [gps_ephem] = alm2geph(alm)
%
% This function converts from an almanac format to a GPS ephemeris format. 
%
% Input:
%   alm        - Yuma almanac format (nx13), with columns of
%                 [ID, health, e, t0, i, asc_rate, sqrt(a), long. of asc_node
%                 at GPS week epoch, perigee, M, af0, af1, week]  
%                 units are seconds, meters, and radians
% Output:
%    gps_ephem - ephemeris matrix for all satellites (nx24), with columns of
%                 [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                 i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                 GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                 Ephemeris parameters are from ICD-GPS-200 with the 
%                 exception of perigee_rate. The gps_ephem will be 
%                 filled with inf values for any almanac element set that 
%                 is out bounds.  
%
% See also: READYUMA, KEP2GEPH, PROPGEPH

% Written by: Jimmy LaMance 10/18/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

% WGS-84 constants
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gps_ephem=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on ALM2GEPH for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'ALM2GEPH';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'alm';
  estruct.variable(1).req_dim = [901 13];
  estruct.variable(1).var = alm;
  
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

% allocate the gps_ephem matrix
gps_ephem = ones(size(alm,1),22) * inf;
                                                      
gps_ephem(:,1) = alm(:,1);                % SV numbers
gps_ephem(:,2) = alm(:,10);               % Mean anomaly
gps_ephem(:,3) = zeros(size(alm(:,1)));   % mean motion correction
gps_ephem(:,4) = alm(:,3);                % eccentricity
gps_ephem(:,5) = alm(:,7);                % sqrt semi-major axis
gps_ephem(:,6) = alm(:,8);                % acsending node
gps_ephem(:,7) = alm(:,5);                % inclination
gps_ephem(:,8) = alm(:,9);                % perigee
gps_ephem(:,9) = alm(:,6);                % ra_rate
gps_ephem(:,10) = zeros(size(alm(:,1)));  % i_rate
gps_ephem(:,11) = zeros(size(alm(:,1)));  % Cuc
gps_ephem(:,12) = zeros(size(alm(:,1)));  % Cus
gps_ephem(:,13) = zeros(size(alm(:,1)));  % Crc
gps_ephem(:,14) = zeros(size(alm(:,1)));  % Crs
gps_ephem(:,15) = zeros(size(alm(:,1)));  % Cic
gps_ephem(:,16) = zeros(size(alm(:,1)));  % Cis
gps_ephem(:,17) = alm(:,4);               % toe (time of ephemeris)
gps_ephem(:,18) = zeros(size(alm(:,1)));  % IODE
gps_ephem(:,19) = alm(:,13);              % GPS week
gps_ephem(:,20) = alm(:,4);               % Toc (time of clock)
gps_ephem(:,21) = alm(:,11);              % Af0 (clock term)
gps_ephem(:,22) = alm(:,12);              % Af1 (clock term)
gps_ephem(:,23) = zeros(size(alm(:,1)));  % Af2 (clock term)
gps_ephem(:,24) = zeros(size(alm(:,1)));  % perigee rate (J2)

%%%%% END ALGORITHM CODE %%%%%

% end ALM2GEPH





