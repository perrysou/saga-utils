function writemov(M, file_name);

% writemov(M, file_name);
%
% Play a Matlab movie and write the movie to an MPEG movie file while  
% retaining the original figure window size and color map properties.
%
%   Input:
%        M         - Matlab Movie matrix
%        file_name - file name to write the MPEG Movie (string)
%   Output:
%        none
%
% See also PLAYMOV, WRITEMPG 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: WRITEMPG

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on WRITEMOV for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% get the figure size from the movie matrix
WidthHeight = M(1:2,1)';                   

% generate a new figure
NewFigure = figure;    

% get the position of the new figure in pixels
set(NewFigure,'units','pixels');
NewFigurePos = get(NewFigure,'position');

% set the width and height of the new figure so it is the same as used
% to create the move
set(NewFigure,'position',[NewFigurePos(1:2) WidthHeight]);

% play the movie and get the colormap
movie(NewFigure,M(:,1))
cmap = colormap;

% close the movie window
close(NewFigure);

% write the movie to an MPEG file with the correct color map
writempg(M, cmap, file_name);

%%%%% END ALGORITHM CODE %%%%%

% end of WRITEMOV
