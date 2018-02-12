%
% gpsconst.m
% 
% assigns values for commonly used GPS constants
%

% References: 	Parkinson, et. al., GPS Theory and Applications, V. 1,
%		AIAA, 1996.
% AJHansen 12 May 1997			Initial coding

c            = 299792458.0;          % velocity of light,  m/sec

f_1          = 1575.42e6;            % L1 frequency, Hz
f_2          = 1227.60e6;            % L2 frequency, Hz
 
lambda_1     = c/f_1;                % L1 wavelength, m
lambda_2     = c/f_2;                % L2 wavelength, m
 
R_sv         = 26561750;             % SV orbit semimajor axis, m

mu_e         = 3.986005e14;          % Earth's grav. parameter (m^3/s^2)
w_e          = 7292115.1467e-11;     % Earth's angular velocity (rad/s)
a_e          = 6378137;              % Earth's semimajor axis, m 
b_e          = 6356752.314;          % Earth's semiminor axis, m 
R_e          = a_e;                  % Earth's approximate radius, m

H_iono       = 350000;               % altitude of the ionospher, m
R_iono       = R_e + H_iono;         % Iono's approximate radius, m
gamma        = (f_1/f_2)^2;          % ionospheric constant for L1/L2
kTEC         = f_1^2*f_2^2/(f_1^2+f_2^2)/40.3;

leap_seconds = 12;
sec_per_day  = 24*3600;
sec_per_wk   = 7*sec_per_day;
TECU2l1m     = 0.163;
