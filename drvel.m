function [v, theta, ar, psi, vc] = drvel(a, h, b, f, g)
% given a h b f g, compute drift velocity, axial ratio and orientation
%   0<=alpha<=1

v = sqrt((g .* h - f .* b).^2+(f .* h - g .* a).^2) ./ (a .* b - h.^2);
theta = atan2(f.*h-g.*a, g.*h-f.*b);

alpha = 2 * sqrt(a.*b-h.^2) ./ (a + b);
ar = sqrt((1 + sqrt(1-alpha.^2))./(1 - sqrt(1-alpha.^2)));
psi = atan2(2*h, a-b) ./ 2 * 180 / pi;

Vc0 = sqrt((a .* b - h.^2)./(a .* g.^2 - 2 * f .* g .* h + b .* f.^2)-1) .* v;
Vc1 = sqrt((a.^3 + 2 * a .* h.^2 + b .* h.^2)./(g .* h + a .* f).^2-1) .* v;


% Vc2 = sqrt((b^3+2*b*h^2+a*h^2)/(f*h+b*g)^2-1)*V;
% Vc0 = (Vc1+Vc2)/2;

if ar >= 6
    vc = Vc1;
else
    vc = Vc0;
end

vbar = mean(v);
thetabar = mean(theta);
arbar = mean(ar);
psibar = mean(psi);
vcbar = mean(vc);


return;
if Vc >= V
    disp('Warning! Characteristic Velocity is larger than Apparent Velocity');
    Vc = NaN;
elseif imag(Vc) == abs(Vc)
    disp('Warning! Characteristic Velocity is imaginary')
    Vc = NaN;
end


end