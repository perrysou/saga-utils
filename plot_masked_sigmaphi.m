function [spm, s4m] = plot_masked_sigmaphi(DATA, DATA_el, kk, dt, elmask, color, tlim)
sv = find(DATA(:, 4) == kk);
if ~isempty(sv)
    svdata = DATA(sv, :);
    ortw = DATA(sv, 1);
    orts = DATA(sv, 2);
    time_sc = gps2utc([ortw, orts]);
    time_sc = datenum(time_sc);
    sv_el = find(DATA_el(:, 4) == kk & DATA_el(:, 3) > elmask);
    if ~isempty(sv_el)
        ortw_el = DATA_el(sv_el, 1);
        orts_el = DATA_el(sv_el, 2);
        time_el = gps2utc([ortw_el, orts_el]);
        time_el = datenum(time_el);
        
        t_d = find(diff(time_sc) > dt)';
        t_d_el = find(diff(time_el) > dt)';
        seq_d = [time_sc(1); time_sc(t_d); time_sc(t_d+1); time_sc(end)];
        seq_d = sortrows(seq_d);
        seq_d_el = [time_el(1); time_el(t_d_el); time_el(t_d_el+1); time_el(end)];
        seq_d_el = sortrows(seq_d_el);
        trx{1} = seq_d';
        trx{2} = seq_d_el';
        
        spm = [];
        s4m = [];
        t = find_common_times(trx);
        if ~isempty(t)
            for tt = 1:2:size(t, 1)
                svdata_masked = svdata(time_sc <= t(tt+1) & time_sc >= t(tt), :);
                if ~isempty(svdata_masked)
                    ortw_masked = svdata_masked(:, 1);
                    orts_masked = svdata_masked(:, 2);
                    time_e = gps2utc([ortw_masked, orts_masked]);
                    time_e = datenum(time_e);
                    sigmaphi_e = svdata_masked(:, 3);
                    s4_e = svdata_masked(:, 5);
                    spm = [spm; time_e, sigmaphi_e, ones(size(time_e)) * kk];
                    s4m = [s4m; time_e, s4_e, ones(size(time_e)) * kk];
                end
            end
            %             %plot the masked sigma_phi
            %             hold on;
            %             [time_e,sigmaphi_e] = discont_proc(spm(:,1),spm(:,2),dt);
            %             plot(time_e,sigmaphi_e,'LineWidth',1,'Color',color);
            %             if isempty(tlim)
            %                 datetick('x','HH','keeplimits');
            %             else
            %                 xlim(tlim);
            %                 datetick('x','HH','keeplimits');
            %             end
            %             hold on;
        end
    else
        spm = [];
        s4m = [];
    end
else
    spm = [];
    s4m = [];
end
end
