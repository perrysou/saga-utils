function [tmat, combos] = find_common_times_v2(trx, debug_flag)
%computes the intersection of a group of (time) segments

switch nargin
    case 0
        debug_flag = 1;
        %data example
        trx{1} = [1, 2, 7, 9, 11, 15, 19, 20];
        trx{2} = [8, 12, 16, 17, 18, 20];
        trx{3} = [1, 3, 8, 9, 13, 17];
        trx{4} = [2, 3, 6, 11, 12, 14, 18, 19];
        trx{5} = [1, 2];
        trx{6} = [];
    otherwise
        debug_flag = 0;
end

%remove empty segment(s)
trx_old = trx;
DRX = {};
for iii = 1:size(trx, 2)
    if ~isempty(trx{iii})
        DRX = [DRX, trx{iii}];
    end
end
trx = DRX;

%reshape the segments with rows of heads and tails
trx_ = cell(size(trx));
for ii = 1:size(trx, 2)
    trx_{ii} = reshape(trx{ii}, 2, []);
end

if debug_flag
    %for plotting purposes
    for ii = 1:size(trx_, 2)
        for jj = 1:size(trx_{ii}, 2)
            plot(trx_{ii}(:, jj), ones(size(trx_{ii}(:, jj)))*(length(trx_) + 1 - ii), '.-k');
            hold on;
        end
        text(-2, length(trx_)+1-ii, ['rx', num2str(ii)]);
    end
end

%number of receivers, from all to a pair
for tol = length(trx):- 1:2
    combos{length(trx)-tol+1,:} = nchoosek(1:length(trx), tol);
end

tmat = cell(length(combos));
for numrx = 1:length(combos)
    for numc = 1:size(combos{numrx,:}, 1)
        rxnum = combos{numrx,:}(numc, 1:end);
        if size(trx, 2) > 1
            for ind = 1:length(rxnum) - 1
                n = 1;
                if ind == 1
                    t1 = trx_{rxnum(ind)};
                else
                    t1 = t;
                end
                t2 = trx_{rxnum(ind+1)};
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
        tmat{numrx, numc} = t;
        clear t;
    end
end
if debug_flag == 0
    return;
else
    %for plotting purposes
    t_ = [t(1:2:end)'; t(2:2:end)'];
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
end


