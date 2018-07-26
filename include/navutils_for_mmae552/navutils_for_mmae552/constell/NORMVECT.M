function [x_norm,x_mag] = normvect(x)

% [x_norm, x_mag] = normvect(x);
%
% Function to normalize n-dimensional vectors to have length of 1.
%
% Input:
%   x      - vector to be normalized  (n x m)
% Output:
%   x_norm - normalized vector, magnitude = 1, same direction as x (n x m)
%   x_mag  - magnitude of the vector (n x 1)

% Written by: Jimmy LaMance 10/25/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% No error checking.  Inputs can be vector, scalar, positive, or negative.
%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% initialize the output variable to be all zeros
x_norm = zeros(size(x));

sum = zeros(size(x,1),1);

for i = 1:size(x,2)
  sum = sum + x(:,i).^2;
end % for i = 1:size(x,2)

x_mag = sqrt(sum);

% find all of the magnitudes with non-zero length
I_non_zero = find(x_mag ~= 0);                   

if any(I_non_zero)
  for i = 1:size(x,2)
    x_norm(I_non_zero,i) = x(I_non_zero,i) ./ x_mag(I_non_zero);
  end % for i = 1:size(x,2)
end % if any(I_non_zero)

%%%%% END ALGORITHM CODE %%%%%

% end of NORMVECT

