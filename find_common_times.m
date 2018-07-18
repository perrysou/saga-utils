function [t] = find_common_times(trx)
%computes the intersection of a group of (time) segments
% %data example
% rng default;
% trx{1} = sort(randi(20,1,10));
% trx{2} = sort(randi(20,1,10));
% trx{3} = sort(randi(20,1,10));
% trx{4} = sort(randi(20,1,10));

%remove empty segment(s)
trx_old = trx;
DRX = {};
for iii = 1:size(trx, 2)
    if ~isempty(trx{iii})
        DRX = [DRX, trx{iii}];
    end
end
trx = DRX;
trx{:};
%reshape the segments with rows of heads and tails
for ii = 1:size(trx, 2)
    tt = trx{ii};
    tstt = tt(1:2:end);
    tend = tt(2:2:end);
    trx_{ii} = [tstt; tend];
end

% %for plotting purposes
% for ii = 1:size(trx_,2)
%     for jj = 1:size(trx_{ii},2)
%         plot(trx_{ii}(:,jj),ones(size(trx_{ii}(:,jj)))*(length(trx_)+1-ii),'.-k');
%         hold on;
%     end
%     text(-2,length(trx_)+1-ii,['rx' num2str(ii)]);
% end

if size(trx, 2) > 1
    for rxnum = 1:size(trx, 2) - 1
        n = 1;
        if rxnum == 1
            t1 = trx_{rxnum};
        else
            t1 = t;
        end
        t2 = trx_{rxnum+1};
        for jj = 1:size(t1, 2)
            for kk = 1:size(t2, 2)
                if t2(1, kk) < t1(1, jj)
                    if t2(2, kk) >= t1(1, jj) && t2(2, kk) <= t1(2, jj)
                        t(1, n) = t1(1, jj);
                        t(2, n) = t2(2, kk);
                    elseif t2(2, kk) > t1(2, jj)
                        t(1, n) = t1(1, jj);
                        t(2, n) = t1(2, jj);
                    else
                        t(1, n) = NaN;
                        t(2, n) = NaN;
                    end
                elseif t2(1, kk) >= t1(1, jj) && t2(1, kk) <= t1(2, jj)
                    t(1, n) = t2(1, kk);
                    if t2(2, kk) <= t1(2, jj)
                        t(2, n) = t2(2, kk);
                    else
                        t(2, n) = t1(2, jj);
                    end
                else
                    t(1, n) = NaN;
                    t(2, n) = NaN;
                end
                n = n + 1;
            end
        end
        t = t(~isnan(t));
        t = reshape(t, 2, []);
    end
    t = t(~isnan(t));
elseif size(trx, 2) == 1 && size(trx_old, 2) == 1
    t = trx{1}';
else
    t = [];
end

return;

%for plotting purposes
tstt = t(1:2:end);
tend = t(2:2:end);
t_ = [tstt'; tend'];
for tt = 1:size(t_, 2)
    plot(t_(:, tt), zeros(size(t_(:, tt))), '.-r');
    hold on;
    plot([t_(1, tt), t_(1, tt)], [0, length(trx_)], '--k');
    plot([t_(2, tt), t_(2, tt)], [0, length(trx_)], '--k');
end
text(-2, 0, 'overlaps');
% plot(rand(20,1)+5,'k');
% axis([-4 20 0 8]);
set(gca, 'yticklabel', []);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
% set(gca,'Visible','off');
% op_path = '/data1/home/ysu27/Dropbox/';
% op_path = '/home/yang/Dropbox/research-misc/ION_2015/';
% saveas(gcf,[op_path,'findcommon.eps'],'epsc2');
close;
