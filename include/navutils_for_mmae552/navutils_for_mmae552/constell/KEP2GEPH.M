function [gps_ephem] = kep2geph(kep_elems,j2_flag)

% [gps_ephem] = kep2geph(kep_elems,j2_flag);
%
% This function converts from a 6-element Keplerian element set
% to a GPS ephemeris format (24 column format). 
%
% Input:
%   kep_elems - Keplerian element set at epoch (nx9), with columns of
%                [sv_num,a,e,i,long. of asc_node,perigee,M,GPS_week,GPS_sec]
%                units are seconds, meters, and radians (mod 2*pi)
%   j2_flag   - flag to indicate the application of J2 effects on
%                ascending node, mean motion, and argument of
%                perigee (optional)
%                0 - no J2 effect
%                1 - apply J2 terms (default)
% Output:
%   gps_ephem - ephemeris matirx for all satellites (nx24), with columns of
%                [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                Ephemeris parameters are from ICD-GPS-200 with the 
%                exception of perigee_rate. The gps_ephem will be 
%                filled with inf values for any kep_elems element set that 
%                is out bounds.  
%
% See also CONSTELL, PROPGEPH 

% See "Global Positioning System: Theory and Applications
% Vol. 1", Spilker et al. or ICD-GPS-200 for details on the meaning
% of the ephemeris variables.

% Written by: Jimmy LaMance 10/18/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

% WGS-84 constants
RADIUS_EARTH = 6378137;       % meters WGS-84 value
MU_EARTH = 3.986005e14;       % m^3/s^2 WGS-84 including the Earth's atmosphere
EARTH_RATE = 7.2921151467e-5; % WGS-84 value in rad/s 
J2 = .00108263;               % IAU 1981

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gps_ephem=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on KEP2GEPH for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% set the j2_flag = 1 (default) if only one input is given
if nargin == 1
  j2_flag = 1;
end % if nargin == 1  

estruct.func_name = 'KEP2GEPH';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'kep_elems';
estruct.variable(1).req_dim = [901 9];
estruct.variable(1).var = kep_elems;
  
estruct.variable(2).name = 'j2_flag';
estruct.variable(2).req_dim = [1 1];
estruct.variable(2).var = j2_flag;
                                
estruct.variable(3).name = 'kep_elems';
estruct.variable(3).req_dim = [901 2; 1 2];
estruct.variable(3).var = kep_elems(:,2:3);
estruct.variable(3).type = 'KEPLER_ORBIT';

% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%% 

% allocate the gps_ephem matrix
gps_ephem = ones(size(kep_elems,1),22) * inf;

% fill in the ephemeris matrix
gps_ephem(:,1) = kep_elems(:,1);            % SV numbers
gps_ephem(:,2) = kep_elems(:,7);            % Mean anomaly 

gps_ephem(:,4) = kep_elems(:,3);            % eccentricity
gps_ephem(:,5) = sqrt(kep_elems(:,2));      % sqrt of the semi-major axis
gps_ephem(:,6) = kep_elems(:,5);            % long. of ascending node
% correct node to comform to the GPS format which has an Earth rate
% correction with toe (see ICD-GPS-200 equations)
%       node         =       node           +  EARTH_RATE *       toe
gps_ephem(:,6) = gps_ephem(:,6) +  EARTH_RATE * kep_elems(:,9);

% check that the new node value is mod 2 * pi
I_node_big = find(gps_ephem(:,6) > 2 * pi);
if any(I_node_big)
  gps_ephem(I_node_big,6) = ...
            rem(gps_ephem(I_node_big,6), 2 * pi);
end % if any(I_node_big)
           
gps_ephem(:,7) = kep_elems(:,4);                    % inclination
gps_ephem(:,8) = kep_elems(:,6);                    % perigee

a = kep_elems(:,2);      % temporary assignment for ease of writing eqs
mean_motion = sqrt(MU_EARTH ./ a.^3);  % sqrt(MU_EARTH/a^3);

% compute right ascension, mean motion, and perigee rates from J2 
if j2_flag == 0
  gps_ephem(:,9) = zeros(size(kep_elems(:,1)));  % ra_rate
  gps_ephem(:,3) = zeros(size(kep_elems(:,1)));  % mean motion corr
  gps_ephem(:,24) = zeros(size(kep_elems(:,1))); % perigee corr
else  
  % J2 contribution terms from p. 64 Hofmann-Wellenhof
  % compute J2 ascending node correction term (m-dot)
  gps_ephem(:,9) = -J2 .* 1.5 .* mean_motion .* RADIUS_EARTH^2 .* ...
                   cos(gps_ephem(:,7)) ./ ...
                   (a.^2 .* (1 - gps_ephem(:,4).^2).^2);
  % compute J2 mean anomaly correction term (m-dot)
  gps_ephem(:,3) = J2 .* .75 .* mean_motion .* RADIUS_EARTH^2 .* ...
                 (3 .* cos(gps_ephem(:,7)).^2 - 1) ./ ...
                 (a.^2 .* sqrt(1 - gps_ephem(:,4).^2).^3);
  % compute J2 argument of perigee correction term (omega-dot)
  gps_ephem(:,24) =  J2 .* .75 .* mean_motion .* RADIUS_EARTH^2 .* ...
                 (5 .* cos(gps_ephem(:,7)).^2 - 1) ./ ...
                 (a.^2 .* (1 - gps_ephem(:,4).^2).^2);
end % if j2_flag == 0

% zero out parameters that are not computed
gps_ephem(:,10) = zeros(size(kep_elems(:,1)));  % i_rate
gps_ephem(:,11) = zeros(size(kep_elems(:,1)));  % Cuc
gps_ephem(:,12) = zeros(size(kep_elems(:,1)));  % Cus
gps_ephem(:,13) = zeros(size(kep_elems(:,1)));  % Crc
gps_ephem(:,14) = zeros(size(kep_elems(:,1)));  % Crs
gps_ephem(:,15) = zeros(size(kep_elems(:,1)));  % Cic
gps_ephem(:,16) = zeros(size(kep_elems(:,1)));  % Cis
gps_ephem(:,18) = zeros(size(kep_elems(:,1)));  % IODE
gps_ephem(:,20) = zeros(size(kep_elems(:,1)));  % Toc (time of clock)
gps_ephem(:,21) = zeros(size(kep_elems(:,1)));  % Af0 (clock term)
gps_ephem(:,22) = zeros(size(kep_elems(:,1)));  % Af1 (clock term)
gps_ephem(:,23) = zeros(size(kep_elems(:,1)));  % Af2 (clock term)

% assign the clock (timing) terms
gps_ephem(:,17) = kep_elems(:,9);               % toe
gps_ephem(:,19) = kep_elems(:,8);               % GPS week

%%%%% END ALGORITHM CODE %%%%%

% end KEP2GEPH




