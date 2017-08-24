function [newtslist, newtelist] = dividet_v3(t, dtau, tolw)
% load tdata.mat
tslist = t(1:2:end);
telist = t(2:2:end);
%length of each continuous segment
tl = telist - tslist;

newtslist = [];
newtelist = [];
for ti = 1:length(tl)
    if tl(ti) >= dtau
        res = mod(tl(ti), dtau);
        quo = (tl(ti) - res) / dtau;
        if res <= tolw
            ddtau = tl(ti) / quo;
            quo1 = quo;
        else
            ddtau = tl(ti) / (quo + 1);
            quo1 = quo + 1;
        end
        %         for numquo = 1:quo1
        %             newtslist = [newtslist; tslist(ti) + (numquo - 1) * ddtau];
        %             newtelist = [newtelist; tslist(ti) + numquo * ddtau];
        %         end
        newtslist = [newtslist; (tslist(ti):tolw:telist(ti) - dtau)'];
        newtelist = [newtelist; (tslist(ti) + dtau:tolw:telist(ti))'];
        if telist(ti) - dtau - newtslist(end) >= tolw / 2
            newtslist = [newtslist; telist(ti) - dtau];
            newtelist = [newtelist; telist(ti)];
        else
            newtelist(end) = telist(ti);
        end
    elseif tl(ti) >= 10
        newtslist = [newtslist; tslist(ti)];
        newtelist = [newtelist; telist(ti)];
    end
end
[newtslist, newtelist, newtelist - newtslist]
end

