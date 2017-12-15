function [tecveq, of_] = stec2vtec(tecl, tecp, el, prns)
bsat_all = getDCB();
bsat = bsat_all(prns)';
brx = -17;
RE = 6.3781 * 10^6;
h = 350 * 1000;
B_rs = sum(tecp - tecl.* sin(el).^2, 'omitnan') ./ sum(sin(el).^2, 'omitnan');
B_rs(B_rs==0) = NaN;
tec = tecl + repmat(B_rs, size(tecl, 1), 1);
tecs = tec - repmat(bsat, size(tec, 1), 1) - brx;

of = cos(asin(RE / (RE + h) * sin(pi / 2 + el)));
of_ = cos(asin(RE / (RE + h) * cos(el)));

tecveq = tecs .* of_;
end