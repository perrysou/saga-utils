function fig_handle = playmov(M, num_plays, frames_per_second, color_map);

% fig_handle = playmov(M, num_plays, frames_per_second, color_map);
%
% Plays Matlab movies in a new figure window while retaining
% the original window size and optionally the original color map. 
%
% Input:
%   M                 - Matlab Movie matrix
%   num_plays         - number of time to play through the movie (optional) 
%                        (1x1), default = 1
%   frames_per_second - number of frames per second to play the movie (optional)
%                        (1x1), default = 2
%   color_map         - color map used for this movie (optional) 
%                        default = Matlab default color map
% Output:
%   fig_handle        - figure handle to the resulting figure
%
% See also WRITEMOV, WRITEMPG 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
fig_handle=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(1,4,nargin);
if ~isempty(msg)
  fprintf('%s  See help on PLAYMOV for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

% verify that num_plays is 1x1 
if nargin < 2
  num_plays = 1;      % set the default value
end %  if nargin > 1

% verify that frames_per_second is 1x1 
if nargin < 3
  frames_per_second = 2;      % set the default value
end %  if nargin > 1   

% sanity checking on the number of plays and the frames per second
if num_plays < 1
  num_plays = 1;
  fprintf('The input number of times to play the movie was less than 1.\n');
  fprintf('The num_plays variable has been to the default value of 1.\n');           
end % if num_plays < 1  

if frames_per_second < 0
  frames_per_second = 2;
  fprintf('The input number of frames per second was less than 0.\n');
  fprintf('The frames_per_second variable has been to the default value of 2.\n');           
end % if frames_per_second < 1  

% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0                        
  estruct.func_name = 'PLAYMOV';

  % Develop the error checking structure with required dimension, matching
  % dimension flags, and input dimensions.
  estruct.variable(1).name = 'num_plays';
  estruct.variable(1).req_dim = [1 1];
  estruct.variable(1).var = num_plays;
  
  estruct.variable(2).name = 'frames_per_second';
  estruct.variable(2).req_dim = [1 1];
  estruct.variable(2).var = frames_per_second;
  
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

% get the figure size from the movie matrix
WidthHeight = M(1:2,1)';                   

% generate a new figure at the bottom left of the screen
NewFigure = figure('units','pixels','position',[50 50 400 400]);    

% get the position of the new figure in pixels
NewFigurePos = get(NewFigure,'position');

% set the width and height of the new figure so it is the same as used
% to create the move
set(NewFigure,'position',[NewFigurePos(1:2) WidthHeight]);

% if a color map was supplied, apply it to the current figure 
if nargin > 3
  colormap(color_map)
end % if nargin > 3

% actually play the movie using the Matlab movie function
movie(NewFigure,M,num_plays,frames_per_second)

%%%%% END ALGORITHM CODE %%%%%

% end of PLAYMOV
