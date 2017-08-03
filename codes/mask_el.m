function [SPM] = mask_el(outdir,signal,rcvr_name,mmddyy_beg, HOUR_init, HOUR_end, DATA, DATA_el)
%%
KK = unique(DATA(:,4));
SPM = [];
switch rcvr_name
    case 'grid108'
        color = 'b';
    case 'grid112'
        color = 'g';
    case 'grid160'
        color = 'k';
end
        
for ki = 1:length(KK)
    kk = KK(ki);
    sv = find(DATA(:,4)==kk);
    svdata = DATA(sv,:);
    ortw = DATA(sv,1);
    orts = DATA(sv,2); 
    time_sc = gps2utc([ortw orts],0);
    time_sc = (time_sc(:,4)) + time_sc(:,5)/60 + time_sc(:,6)/3600;
    
    %remove data from the next day(temp)
    [ ind ] = same_day_check( time_sc );
    time_sc = time_sc(ind);
     
    
    %plot el
    sv_el = find(DATA_el(:,4)==kk);
    ortw_el = DATA_el(sv_el,1);  
    orts_el = DATA_el(sv_el,2);
    el = DATA_el(sv_el,3);
    time_el = gps2utc([ortw_el orts_el],0);
    time_el = (time_el(:,4)) + time_el(:,5)/60 + time_el(:,6)/3600;
    dt = 0.5;
    if size(time_el,1)>4
        [time_el_e,el_e] = discont_proc(time_el,el,dt);
        if kk<=16
            fig1 = subplot(8,2,kk);   
        elseif kk == 17
            figure;
            subplot(8,2,kk-16);
        else
            fig2 = subplot(8,2,kk-16);
        end

        %remove data from the next day(temp)
        [ ind ] = same_day_check( time_el_e );
        time_el_e = time_el_e(ind);
        el_e = el_e(ind);
        
        set(gca,'FontSize',8,'YTick',1/90*[30 90],...,
        'YTickLabel',{'30/0.066','90/0.2'},'YGrid','on');
        plot(time_el_e,el_e/90,'LineWidth',1,'Color',color);
        set(gca,'XTick',0:6:24,'XGrid','on');
        xlim([0 24]);    
        ylim([0 1.1]);
        xlabel(['PRN',num2str(kk)]);  
    end
    
    %plot masked sigmaphi
    sv_el = find(DATA_el(:,4)==kk & DATA_el(:,3)>30);
    ortw_el = DATA_el(sv_el,1);  
    orts_el = DATA_el(sv_el,2);
    spm = [];
    
    if isempty(orts_el) || isempty(ortw_el)== 1
        disp(['PRN ',num2str(kk),' is at low elevation(<30 degrees).']);
    else
        time_el = gps2utc([ortw_el orts_el],0);
        time_el = (time_el(:,4)) + time_el(:,5)/60 + time_el(:,6)/3600;
        if size(time_sc,1)>4
            dt = mean(diff(time_sc));
            t_d = find(diff(time_sc)>0.5)';
            
            %remove data from the next day(temp)
            [ ind ] = same_day_check( time_el );
            time_el = time_el(ind);
            
            t_d_el = find(diff(time_el)>0.5)';                   
            seq_d = [time_sc(1); time_sc(t_d); time_sc(t_d+1); time_sc(end)];
            seq_d = sortrows(seq_d);
            seq_d_el = [time_el(1); time_el(t_d_el); time_el(t_d_el+1); time_el(end)];
            seq_d_el = sortrows(seq_d_el);
            trx{1} = seq_d';
            trx{2} = seq_d_el';
            t = find_common_times(trx);
            for tt = 1:2:size(t,1)
                svdata_masked = svdata(time_sc<=t(tt+1) & time_sc>=t(tt),:);
                ortw_masked = svdata_masked(:,1);orts_masked = svdata_masked(:,2);
                time_e = gps2utc([ortw_masked orts_masked],0);
                time_e = (time_e(:,4)) + time_e(:,5)/60 + time_e(:,6)/3600;
                sigmaphi_e = svdata_masked(:,3);
                SPM = [SPM;time_e sigmaphi_e ones(size(time_e))*kk];
                spm = [spm;time_e sigmaphi_e];
            end
            hold on;
            [time_sc_e,sigmaphi_e] = discont_proc(spm(:,1),spm(:,2),dt);
            plot(time_sc_e,sigmaphi_e*5,'LineWidth',1,'Color',color);           
        end    
    if kk == KK(1)||kk == KK(17)
%     timestr = strcat({'Time [hours] after '},num2str(HOUR_init),{':00 UT on (mm/dd/yy) '},mmddyy_beg);
%     titlestr = strcat({'Low rate masked S_4 with elevation>30 degrees for receiver: '},rcvr_name);
%     text(6,-26,0,[timestr titlestr],'Fontsize',8);
    end
    end   
end

%sort sigmaphi data by timestamps
SPM = sortrows(SPM,1);
SPM(:,2) = SPM(:,2)*2*pi;
fig_name = strcat(signal,'_LowRateSigmaPhiWithElevationMask_',rcvr_name);
fig_1 = strcat(outdir,fig_name,'_fig1','.png');
fig_2 = strcat(outdir,fig_name,'_fig2','.png');
saveas(fig1,fig_1);
saveas(fig2,fig_2);
end

