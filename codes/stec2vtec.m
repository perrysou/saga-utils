function [tecveq, of_] = stec2vtec(tecl, tecp, el, prns)
bsat_all = getDCB();
bsat = bsat_all(prns)';
brx = 0;
RE = 6.3781 * 10^6;
h = 350 * 1000;
B_rs = sum(tecp - tecl .* sin(el), 'omitnan') ./ sum(sin(el), 'omitnan');
B_rs(B_rs==0) = NaN;
tec = tecl + B_rs;
tecs = tec - bsat - brx;

of = cos(asin(RE / (RE + h) * sin(pi / 2 + el)));
of_ = cos(asin(RE / (RE + h) * cos(el)));

tecveq = tecs .* of_;
end