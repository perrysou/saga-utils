function [newtslist, newtelist] = dividet_v1(t, dtau, tolw)
tslist = t(1:2:end);
telist = t(2:2:end);

%length of each continuous segment
tl = telist - tslist;

newtslist = [];
newtelist = [];
for ti = 1:length(tl)
    if tl(ti) >= 30
        res = mod(tl(ti), dtau);
        quo = (tl(ti) - res) / dtau;
        if res <= tolw
            ddtau = tl(ti) / quo;
            quo1 = quo;
        else
            ddtau = tl(ti) / (quo + 1);
            quo1 = quo + 1;
        end
        for numquo = 1:quo1
            newtslist = [newtslist; tslist(ti) + (numquo - 1) * ddtau];
            newtelist = [newtelist; tslist(ti) + numquo * ddtau];
        end
    end
end
end