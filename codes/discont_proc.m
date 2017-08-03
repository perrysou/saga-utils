function [time_e,data_e] = discont_proc(time,data,dt)
%Identify discontinuities in time/data
% time = [1,2,4,5,6,8,13]';
% data = rand(7,1);
% dt = 2;
t_d = find(diff(time)>dt)';
seq_d = [1 t_d length(time)];
time_e = [];data_e = [];
        if isempty(t_d)
            time_e = time;
            data_e = data;
        else
            for ind_d = 1:1:length(t_d)
                time_d = (time(seq_d(ind_d+1))+dt:dt:time(seq_d(ind_d+1)+1))';
                time_d = time_d*NaN';
                data_d = ones(size(time_d))*NaN;
                if ind_d~=1
                   time_e = [time_e;time(seq_d(ind_d)+1:seq_d(ind_d+1));time_d];
                    data_e = [data_e;data(seq_d(ind_d)+1:seq_d(ind_d+1));data_d];
                else
                    time_e = [time_e;time(seq_d(ind_d):seq_d(ind_d+1));time_d];
                    data_e = [data_e;data(seq_d(ind_d):seq_d(ind_d+1));data_d];
                end
            end
            time_e = [time_e;time(seq_d(end-1)+1:seq_d(end))];
            data_e = [data_e;data(seq_d(end-1)+1:seq_d(end))];
        end
end
