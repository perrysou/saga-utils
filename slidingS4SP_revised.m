function [S4out, SPHout] = slidingS4SP_revised(Powin, PHin, N)
% derived from slidingavg
%   OUTPUT_ARRAY = SLIDINGAVG(INPUT_ARRAY, N)
%
%  The function 'slidingavg' implements a one-dimensional filtering,
%  applying a sliding window to a sequence. Such filtering replaces the
%  center value in the window with the average value of all the points
%  within the window. When the sliding window is exceeding the lower or
%  upper boundaries of the input vector INPUT_ARRAY, the average is
%  computed among the available points. Indicating with nx the length of
%  the the input sequence, we note that for values of N larger or equal to
%  2*(nx - 1), each value of the output data array are identical and equal
%  to mean(in).
%
%  *  The input argument INPUT_ARRAY is the numerical data array to be
%  processed. *  The input argument N  is the number of neighboring data
%  points to average over for each point of IN.
%
%  *  The output argument OUTPUT_ARRAY is the output data array.
%
%     � 2002 - Michele Giugliano, PhD and Maura Arsiero
%     (Bern, Friday July 5th, 2002 - 21:10)
%    (http://www.giugliano.info) (bug-reports to michele@giugliano.info)
%
% Two simple examples with second- and third-order filters are
%  slidingavg([4 3 5 2 8 9 1],2)
%  ans =
%   3.5000  4.0000  3.3333  5.0000  6.3333  6.0000  5.0000
%
%  slidingavg([4 3 5 2 8 9 1],3)
%  ans =
%   3.5000  4.0000  3.3333  5.0000  6.3333  6.0000  5.0000
%

if (isempty(Powin)) || (isempty(PHin)) || (N <= 0)
    % If the input array is empty or N is non-positive,% an error is reported
    % to the standard output and the  % execution of the routine is stopped.
    disp(sprintf('SlidingAvg: (Error) empty input data or N null.'));
    return;
    
end % if

%find amplitude
amp = sqrt(Powin);

% if (N==1) If the number of neighbouring points over which the sliding
%  S4out = Powin; SPHout = PHin;
%     out = in;
% average will be performed is '1', then no average actually occur and
%  return;
% OUTPUT_ARRAY will be the copy of INPUT_ARRAY and the execution of the
% routine end % if is stopped.

nx = length(amp); %=length(PHin)
% The length of the input data structure is acquired to later evaluate the
% 'mean' over the appropriate boundaries.


% If the number of neighbouring points over which the sliding average
% will be performed is large enough, then the average actually covers
% all the points
% of INPUT_ARRAY, for each index of OUTPUT_ARRAY and some CPU time can be
% gained by such an approach.
% The execution of the routine is stopped.
if (N >= (2 * (nx - 1)))
    S4out = sqrt((mean(Powin) - mean(amp).^2)/mean(amp).^2) * ones(size(Powin));
    SPHout = std(PHin) * ones(size(PHin)); %sigma_phi
    return;
end % if

S4out = zeros(size(Powin));
% In all the other situations, the initialization of the output data
% structure is performed.
SPHout = zeros(size(PHin));

% When N is even, then we proceed in taking the half of it:% m = N
% /  2.% Otherwise (N >= 3, N odd), N-1 is even ( N-1 >= 2) and we
% proceed taking the half of it:
% m = (N-1) /  2.
if rem(N, 2) ~= 1
    m = N / 2;
else
    m = (N - 1) / 2;
end % if

for i = 1 + m:nx - m,
    % For each element (i-th) contained in the input numerical array, a
    % check must be performed:
    if ((i - m) < 1) && ((i + m) <= nx)
        % If not enough points are available on the left of the i-th element
        % then we proceed to evaluate the mean from the first element to the
        % (i + m)-th.
        S4out(i) = sqrt((mean(Powin(1:i+m)) - mean(amp(1:i+m)).^2)/ ...
            mean(amp(1:i+m)).^2);
        SPHout(i) = std(PHin(1:i+m));
        %   out(i) = mean(in(1:i+m));
    elseif ((i - m) >= 1) && ((i + m) <= nx)
        % If enough points are available on the left and on the right of the
        % i-th element then we proceed to evaluate the mean on 2*m
        % elements centered on the i-th position.
        S4out(i) = sqrt((mean(Powin(i-m:i+m)) - mean(amp(i-m:i+m)).^2)/ ...
            mean(amp(i-m:i+m)).^2);
        SPHout(i) = std(PHin(i-m:i+m));
        %   out(i) = mean(in(i-m:i+m));
    elseif ((i - m) >= 1) && ((i + m) > nx)
        % If not enough points are available on the rigth of the i-th element
        % then we proceed to evaluate the mean from the element (i - m)-th to
        % the last one.
        S4out(i) = sqrt((mean(Powin(i-m:nx)) - mean(amp(i-m:nx)).^2)/ ...
            mean(amp(i-m:nx)).^2);
        SPHout(i) = std(PHin(i-m:nx));
        %   out(i) = mean(in(i-m:nx));
    elseif ((i - m) < 1) && ((i + m) > nx)
        % If not enough points are available on the left and on the rigth of
        % the i-th element then we proceed to evaluate the mean from the first
        % element to the last.
        S4out(i) = sqrt((mean(Powin(1:nx)) - mean(amp(1:nx)).^2)/ ...
            mean(amp(1:nx)).^2);
        SPHout(i) = std(PHin(1:nx));
        %   out(i) = mean(in(1:nx));
    end % if
end % for i
