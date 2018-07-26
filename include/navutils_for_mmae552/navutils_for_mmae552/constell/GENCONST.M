function [constell_elems] = genconst(constell_define)

% [constell_elems] = genconst(constell_define)
%
% This function converts information about satellite constellations
% into the corresponding 6-element Keplerian element sets
%
%   Input:
%       constell_define - Array with each row defining a constellation
%                         [nsats - # of satellites in constellation,
%                          nplanes - # of orbit planes,
%                          a - semimajor axis (meters),
%                          e - eccentricity,
%                          i - inclination (rad),
%                          ASC_i - longitude of ascending node of 1st plane 
%                                  (rad) (Optional),
%                          w_i   - argument of perigee of 1st plane (rad) 
%                                  (Optional)
%                          M_i   - mean anomaly of 1st satellite in 1st plane 
%                                  (rad) (Optional)
%                          delta_ASC - longitude of ascending node spacing between
%                                      adjacent planes (rad) (Optional)
%                          delta_w - argument of perigee spacing between adjacent 
%                                    planes (rad) (Optional)
%                          delta_M - delta mean anomaly between satellites in
%                                    adjacent planes (rad) (Optional)] 
%       Note: If only 5 parameters are input, then ASC_i=0, w_i=0, M_i=0,
%             orbit planes will be evenly spaced, all sv's will have same w,
%             and sv's in adjacent planes will have same M.
%                          
%   Output:
%       constell_elems - Constellation ephemeris [sv, a, e, i, long. of asc node, w, M]
%                        (meters, radians)
%       Note: sv's start counting at 101

% Written by: Maria Evans 12/19/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
constell_elems=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,1,nargin);
if ~isempty(msg)
  fprintf('%s  See help on GENCONST for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'GENCONST';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'constell_define';
  estruct.variable(1).req_dim = [901 5; 901 6; 901 7; 901 8; 901 9; 901 10; ...
                                 901 11;];
  estruct.variable(1).var = constell_define;
  
  estruct.variable(2).name = 'inclination';
  estruct.variable(2).var = constell_define(:,5);
  estruct.variable(2).type = 'ANGLE_RAD';

  var_num = 3;
  if size(constell_define,2) >= 6,
    estruct.variable(var_num).name = 'ASC_i';
    estruct.variable(var_num).var = constell_define(:,6);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  if size(constell_define,2) >= 7,
    estruct.variable(var_num).name = 'w_i';
    estruct.variable(var_num).var = constell_define(:,7);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  if size(constell_define,2) >= 8,
    estruct.variable(var_num).name = 'M_i';
    estruct.variable(var_num).var = constell_define(:,8);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  if size(constell_define,2) >= 9,
    estruct.variable(var_num).name = 'delta_ASC';
    estruct.variable(var_num).var = constell_define(:,9);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  if size(constell_define,2) >= 10,
    estruct.variable(var_num).name = 'delta_w';
    estruct.variable(var_num).var = constell_define(:,10);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  if size(constell_define,2) == 11,
    estruct.variable(var_num).name = 'delta_M';
    estruct.variable(var_num).var = constell_define(:,11);
    estruct.variable(var_num).type = 'ANGLE_RAD';
    var_num = var_num + 1;
  end;

  % Call the error checking function
  stop_flag = err_chk(estruct);
  
  if stop_flag == 1           
    fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
             estruct.func_name);
    return
  end % if stop_flag == 1
end % if matlab_version >= 5.0 

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% fill in the constellation ephemeris matrix
[n_constell,params]=size(constell_define);

ASC_i(1:n_constell)=0;
w_i(1:n_constell)=0;
M_i(1:n_constell)=0;
if params<9,
  delta_ASC(1:n_constell)=2*pi/constell_define(:,2); 
else
  delta_ASC(1:n_constell)=constell_define(:,9);
end;

if params>=6, ASC_i=constell_define(:,6); end;
if params>=7, w_i=constell_define(:,7); end;
if params>=8, M_i=constell_define(:,8); end;
if params>=10,
  delta_w(1:n_constell)=constell_define(:,10);
else,
  delta_w(1:n_constell)=zeros(n_constell,1);
end;
if params==11,
  delta_M(1:n_constell)=constell_define(:,11);
else,
  delta_M(1:n_constell)=zeros(n_constell,1);
end;

constell_elems(:,1) = [101:100+sum(constell_define(:,1))]';
nsats_per_plane = constell_define(:,1)./constell_define(:,2);

sv = 1;
for i=1:n_constell,     % for each constellation
  nsats = constell_define(i,1);
  nplanes = constell_define(i,2);
  satsperplane = nsats/nplanes;

  % Semimajor axis
  constell_elems(sv:sv+nsats-1,2) = ones(nsats,1)*constell_define(i,3); 
  % Eccentricity
  constell_elems(sv:sv+nsats-1,3) = ones(nsats,1)*constell_define(i,4); 
  % Inclination
  constell_elems(sv:sv+nsats-1,4) = ones(nsats,1)*constell_define(i,5);

  index = [[1:nplanes]'*ones(1,satsperplane)]';
  index = reshape(index,size(index,1)*size(index,2),1);

  % Ascending Node
  if delta_ASC(i) ~= 0,
    delta_node = [0:delta_ASC(i):delta_ASC(i)*(nplanes-1)]';
  else,
    delta_node = zeros(nplanes,1);
  end;
  constell_elems(sv:sv+nsats-1,5) = ...
               [ones(nsats,1)*ASC_i(i)]+delta_node(index); 

  % Argument of Perigee
  if delta_w(i) ~= 0,
    delta_argp = [0:delta_w(i):delta_w(i)*(nplanes-1)]';
  else,
    delta_argp = zeros(nplanes,1);
  end;
  constell_elems(sv:sv+nsats-1,6) = ...
               [ones(nsats,1)*w_i(i)]+delta_argp(index); 

  % Mean Anomaly
  if delta_M(i) ~= 0,
    delta_mean = [0:delta_M(i):delta_M(i)*(nplanes-1)]';
  else,
    delta_mean = zeros(nplanes,1);
  end;
  M_sat = 2*pi/satsperplane;
  Mi = M_i(i);
  constell_elems(sv:sv+nsats-1,7) = ...
       reshape((reshape(rem([Mi:M_sat:Mi+(nsats-1)*M_sat],2*pi),...
              satsperplane,nplanes)' + delta_mean(:,ones(1,satsperplane)))', ...
              nsats,1);

  sv = sv + nsats;
end;

constell_elems(:,3:7) = rem(constell_elems(:,3:7),2*pi);

%%%%% END ALGORITHM CODE %%%%%

% end of GENCONST
