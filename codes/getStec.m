function [tecl, tecp] = getStec(l1, l2, p1, p2)

% constants
k = 40.308193;
c = 299792458;
fArr = [1575.42, 1227.60] * 10^6;
lambdaArr = c ./ fArr;
f1 = fArr(1); f2 = fArr(2);
lambda1 = lambdaArr(1); lambda2 = lambdaArr(2);

tecl = 1 / k * (f1 * f2)^2 / (f1^2 - f2^2) * (l1 * lambda1 - l2 * lambda2) / 10^16;

tecp = 1 / k * (f1 * f2)^2 / (f1^2 - f2^2) * (p2 - p1) / 10^16;
end