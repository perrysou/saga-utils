function [tecl_, tecp_] = getStec(l1, l2, p1, p2)
l1(l1 == 0) = NaN;
l2(l2 == 0) = NaN;
p1(p1 == 0) = NaN;
p2(p2 == 0) = NaN;
% constants
k = 40.308193;
c = 299792458;
fArr = [1575.42, 1227.60] * 10^6;
ns2tecu = -1 / ((k / fArr(1)^2 - k / fArr(2)^2) / c * 10^9) / 10^16;
lambdaArr = c ./ fArr;
f1 = fArr(1); f2 = fArr(2);
% lambda1 = lambdaArr(1); lambda2 = lambdaArr(2);

tecl = 1 / k * (f1 * f2)^2 / (f1^2 - f2^2) * (l1 * lambdaArr(1) - l2 * lambdaArr(2)) / 10^16;
tecl_ = (l1 / fArr(1) - l2 / fArr(2)) * 10^9 * ns2tecu;
tecp = 1 / k * (f1 * f2)^2 / (f1^2 - f2^2) * (p2 - p1) / 10^16;
tecp_ = (p2 - p1) / c * 10^9 * ns2tecu;

if isempty(tecl_)
    tecl_ = NaN;
end
if isempty(tecp_)
    tecp_ = NaN;
end
end