function res = angle_truncate(angle, unit)
%%% input: any angle in radians or degrees
%%% return an angle in [-pi, pi] or [-180, 180]
if nargin == 1 || (nargin == 2 && strcmp(unit, 'rad'))
    pie = pi;
elseif nargin == 2 && strcmp(unit, 'deg')
    pie = 180;
else
    ME = MException('myComponent:inputError', ...
        'Incorrect input arguments.');
    throw(ME);
end
res = mod(angle + pie, 2 * pie) - pie;
end