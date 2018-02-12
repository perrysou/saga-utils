function [a_n, im, m] = nearest(a, b)
[m, im] = min(abs(repmat(a, 1, length(b))-repmat(b, 1, length(a))'));
a_n = a(im);
end