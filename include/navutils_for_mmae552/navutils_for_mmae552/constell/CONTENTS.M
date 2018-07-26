% Constellation Toolbox Version 4.01-E   8/5/99
%
% EXAMPLES AND GENERAL SUPPORT FUNCTIONS 
%   contents  -  Display contents of the toolbox.
%   err_chk   -  Input variable error checking.
%   ex_dgps   -  Example script for differential GPS positioning.
%   ex_dop_m  -  Example program for DOP coverage over the Earth with a movie.
%   ex_nmea   -  Example program for reading and plotting NMEA data.
%   grnd2orb  -  Example program demonstrating the visibility and DOPs routines
%                for observers fixed on the Earth
%   normvect  -  Computes the unit vector and the magnitude of a given vector.
%   orb2orb   -  Example function demonstrating the visibility and DOPS routines
%                for orbiting observers.
%
% YUMA ALMANACS, EPHEMERIS FORMATS, AND SATELLITE PROPAGATION
%   alm2geph  -  Converts YUMA almanacs to GPS ephemeris format.
%   find_alm  -  Find the latest almanac files in the Matlab path.
%   genconst  -  Constellation Keplerian element generation.
%   kep2geph  -  Converts Keplerian elements (a,e,i,RA,w,M) to GPS ephemeris format.
%   keplr_eq  -  Solves Kepler's equation using a Newton iteration method.
%   propgeph  -  Computes satellite orbits using GPS ephemeris format.
%   readyuma  -  Function to read in YUMA formatted almanacs.
%   writyuma  -  Create a YUMA formatted almanac file.
%
%
% VISIBILITY
%   grnd2orb  -  Example program demonstrating the visibility and DOPs routines
%                for observers fixed on the Earth
%   los       -  Computes line-of-sight (LOS).
%   num_vis   -  Computes the number of satellites visible.
%   orb2orb   -  Example function demonstrating the visibility and DOPS routines
%                for orbiting observers.
%   passdata  -  Determines pass numbers for data that has passed the masking
%                tests in vis_data. Optionally returns relevant pass information.
%   vis_data  -  Determines visible satellites and other data using mask info.
%
% GRAPHICAL OUTPUT
%   makeplot  -  Function to create the five most used plots, azimuth, elevation, 
%                sky plot, number of visible satellites, and DOPS for an object.
%   orb_anim  -  Satellite trajectory animation.
%   playmov   -  Function for playing a Matlab movie.
%   plotpass  -  Function for plotting data pertaining to passes (az, el, etc.)
%   plotsky   -  Function for plotting azimuth/levation data on a sky plot.
%   writemov  -  Function to write MPG movies to a file.
%   writempg  -  Supporting m function for writemov.
%
% DILUTION OF PRECISION (DOP)
%   dops2err  -  Converts DOPS to position error.
%   ll2dops   -  Computes DOPS given local level inputs.
%   ned2dops  -  Computes DOPS given NED inputs.
%
% POSITION ERROR ANALYSIS AND SIMULATION
%   addambig  -  Adds N integer ambiguities to CPH measurements.
%   clockerr  -  Simulates receiver clock bias and drift.
%   doppler   -  Computes raw Doppler measurments.
%   ionodlay  -  Computes the ionosphere path delay using the GPS model.
%   lsnav     -  Least squares navigation solution using PR and CPH measurements
%   pseudo_r  -  Generates PR measurements.
%   sa_clock  -  Simulates the dither effect of S/A on pseudoranges.
%   sa_eps    -  Simulates the epsilon effect of S/A on pseudoranges.
%   tropdlay  -  Computes the troposphere path delay using the Hopfield model.
%
% DATA PROCESSING
%   nmeanext  -  Function for reading NMEA messages.
%   parsegga  -  Parsing function for NMEA GGA messages.
%   parsegsa  -  Parsing function for NMEA GSA messages.
%   parsegsv  -  Parsing function for NMEA GSV messages.
%   parsnmea  -  Driver parsing function for NMEA messages.
%   readnmea  -  Function to read NMEA data.
%   readyuma  -  Function to read in YUMA formatted almanacs.
%   writyuma  -  Create a YUMA formatted almanac file.
%
% DIFFERENTIAL PROCESSING
%   add_dpr   -  Adds differential correction to PR measurements.
%   diffcorr  -  Computes differential corrections.
%
% COORDINATE TRANSFORMATIONS
%   azel2ned  -  Converts from azimuth and elevation to NED coordinates.
%   body2ned  -  Converts from body coordinates to NED coordinates.
%   ecef2eci  -  Converts from ECEF to ECI coordinates.
%   ecef2ll   -  Converts from ECEF to Local Level (LL) coordinates.
%   ecef2lla  -  Converts from ECEF to lat, lon, and hgt.
%   ecef2ned  -  Converts from ECEF to NED coordinates.
%   eci2ecef  -  Converts from ECI to ECEF coordinates.
%   eci2ll    -  Converts from ECI to Local Level coordinates.
%   ll2ecef   -  Convert from local level to ECEF.
%   ll2eci    -  Convert from local level to ECI.
%   lla2ecef  -  Convert from lat, lon, and hgt to ECEF.
%   ned2azel  -  Converts from NED coordinates to azimuth and elevation.
%   ned2body  -  Converts from NED to body coordinates.
%   ned2ecef  -  Converts from NED to ECEF coordinates.
%   sidereal  -  Computes the Greenwich sidereal time.
%
% TIME TRANSFORMATIONS
%   gps2utc      - Converts from GPS to UTC time formats.
%   gpst2sec     - Convert from [gps_week, gps_sec] to [total_gps_sec].
%   leapsecs.dat - Data file with leap second information.
%   sec2gpst     - Convert from [total_gps_sec] to [gps_week, gps_sec].
%   utc2gps      - Converts from UTC to GPS time formats.
%   utc2leap     - Computes the number of leap seconds at a given UTC time.


 
