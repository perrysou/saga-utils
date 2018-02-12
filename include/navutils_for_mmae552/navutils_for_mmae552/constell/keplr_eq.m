function [E] = keplr_eq(M,e)

% [E] = keplr_eq(M,e);
%
% Function (vectorized) to iteratively solve Kepler's 
% equation for eccentric anomaly.
%
% Input:
%   M - mean anomaly (rad) (nxm)
%   e - eccentricity (dimensionless) (nxm) or (1x1)
% Output:
%   E - eccentric anomaly (rad) (nxm)
%
% See also PROPGEPH 

% Written by: Jimmy Lamance 10/21/96
% Copyright (c) 1998 by Constell, Inc.

% functions called: none

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
E=[];

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,2,nargin);
if ~isempty(msg)
  fprintf('%s  See help on KEPLR_EQ for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

I = find(e < 0);

if any(I)
  fprintf('KEPLR_EQ eccentricites must be greater than 0.\n')
  fprintf('Example of bad eccentricity %f.\n',e(I(1)));
  error('Invalid eccentricities (e) in KEPLR_EQ');
end % if any(I) 

% clear out varible checking
clear I min_e max_e

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% set convergence criterion
tol = 1e-12;      % convergence criterion
max_iter = 10;    % max iteration 

m = 1;            % iteration counter initial value
dE = 1.0;         % set initial value of delta-E to fail convergenge

E = M + (e .* sin(M)) ./ (1 - sin(M + e) + sin(M));    % set initial value of E 

while (any(find(dE > tol)) & m <= max_iter)
  dE = (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
  E = E - dE;         % update E
  m = m + 1;          % increment iteration numbers

  if m > max_iter     % check against maximum iterations allowed
    fprintf('Warning message from KEPLR_EQ ...\n')
    fprintf('Maximum iterations exceeded in solution of Kepler''s Equation.\n')
    fprintf('Results may be invalid.\n\n',e)
    return
  end % if

end % while


%%%%% END ALGORITHM CODE %%%%%

% end of KEPLR_EQ

