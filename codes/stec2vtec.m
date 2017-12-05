function [tecveq] = stec2vtec(tecl, tecp, el)
bsat = 0;
brx = 0;
B_rs = sum(tecp - tecl .* sin(el).^2) / sum(sin(el).^2);
tec = tecl + B_rs;
tecs = tec - bsat - brx;
tecveq = tecs .* cos(asin(0.94902 * cos(el)));
end