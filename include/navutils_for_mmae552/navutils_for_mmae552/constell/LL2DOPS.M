function [dops, t_out, num_sats] = ll2dops(x_ll,t);

% [dops, t_out, num_sats] = ll2dops(x_ll,t);
%
% Function to compute DOP values from Local Level (LL) LOS vectors.  All 
% LOS data that is passed into this function will be used in the DOP 
% computations.   Masking and data editing should be performed prior to 
% passing the Local Level LOS vectors to LL2DOPS.
%
% Input:
%   x_ll     - Visible LOS unit vectors in local level coordinates (nx3)
%   t        - GPS time corresponding to each row in the LOS vectors, 
%               [GPS_week GPS_sec]  (nx2)
%               valid GPS_week values are 1-3640 (years 1980-2050)
%               valid GPS_sec values are 0-604799
% Output:
%   dops     - dilution of precision for each time (kx5) 
%               where k = number of time steps
%               [GDOP, PDOP, HDOP, VDOP, TDOP]
%   t_out    - output time matrix (kx2) [GPS_week GPS_sec]
%   num_sats - number of satellites used in the DOP computation (kx1)
%
%   Note:  If only 3 satellites are found at a given time, the altitude is
%          assumed fixed and only the HDOP and TDOP are filled with values.
%          If fewer than 3 satellites are found at a given time, the default
%          values of inf are returned for all of the DOPS for that time.
%
% See also VIS_DATA, NED2DOPS

% Written by: Jimmy LaMance 10/24/96
% Copyright (c) 1998 by Constell, Inc.

% Reference: Parkinson and Spilker (Blue Book) vol. 1, page 413-414. Modified to 
% operate in NED coordinates.

% functions called: ERR_CHK, NED2DOPS

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% Initialize the output variables
dops=[]; t_out=[]; num_sats=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on LL2DOPS for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.func_name = 'LL2DOPS';

estruct.variable(1).name = 'x_ll';
estruct.variable(1).req_dim = [901 3];
estruct.variable(1).var = x_ll;
  
estruct.variable(2).name = 't';
estruct.variable(2).req_dim = [901 2];
estruct.variable(2).var = t;
estruct.variable(2).type = 'GPS_TIME';
                                
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% this is really the same thing as NED to DOPS, so call NED2DOPS
[dops, t_out, num_sats] = ned2dops(x_ll,t);

%%%%% END ALGORITHM CODE %%%%%

% end of LL2DOPS
