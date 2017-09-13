function [tlim,splim,s4lim] = el_and_masked_sigmaphi(cases_folder,home_dir,rcvr_op,doy,year,signal,spmask,s4mask)
sep = filesep;
%% Plot elevation and elevation-filtered sigmaphi prn by prn for all receivers
SPM = [];S4M = [];
for kk = 1:32
    if kk<=16
    fig1 = subplot(8,2,kk);   
    elseif kk == 17
    figure;
    subplot(8,2,kk-16);
    else
    fig2 = subplot(8,2,kk-16);
    end
    
    legend_count = 0;
    for rr = 1:size(rcvr_op,1)
        rcvr_name = rcvr_op(rr,:);
        [color] = rx_color(rcvr_name);
        scintout = strcat('prn_files_',signal,sep,'scintdata.mat');
        elout = strcat('prn_files_',signal,sep,'txinfodata.mat');
%         if strcmp(case_folder(end-4:end-1),'pfrr')
%             %folder_path for 2013 Poker Flat data 
%             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%         else
%             %folder_path for 2013 Calgary data
%             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);     
%         end
        [~,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
        %load scintdata & txinfodata
        load([op_path,scintout]);load([op_path,elout]);
        if ~isempty(DATA)
        datestt = gps2utc(DATA(1,1:2));
        dateend = gps2utc(DATA(end,1:2));
        elseif ~isempty(DATA_el)
            datestt = gps2utc(DATA_el(1,1:2));
            dateend = gps2utc(DATA_el(end,1:2));
        else
            datestt = [];
            dateend =[];
        end  
%         tlim = find_tlim(datestt,dateend,doy);
        tlim = datenum([datestt;dateend])';
        %plot elevation w/o 30deg mask
        dt = 1/24/3600*600;      
        if ~isempty(DATA_el)
            plot_el(DATA_el,kk,dt,color,tlim);
            legend_count = legend_count + 1;
            RCVRNAME{legend_count,:} = [rcvr_name,' el'];           
        end
        if ~isempty(DATA)&&~isempty(DATA_el)
            elmask = 30;
            [spm{rr},s4m{rr}] = plot_masked_sigmaphi(DATA,DATA_el,kk,dt,elmask,color,tlim);
            legend_count = legend_count + 1;
            RCVRNAME{legend_count,:} = [rcvr_name,' \sigma_\phi'];
        else 
            spm{rr} = [];
            s4m{rr} = [];
        end
    end
    if kk == 32 || kk == 16
%         legend(RCVRNAME);
        ylabel('el [deg], $\sigma_\phi$ [cycle]')
    end
    SPM = [SPM;spm];
    S4M = [S4M;s4m];
end

%% Save elevation-filtered and elevation & sp mask filtered scintillation data
MAXSP = [];
MAXS4 = [];
for rr = 1:size(SPM,2)
    rcvr_name = rcvr_op(rr,:);
    SPM_M = [];
    S4M_M = [];
    for kk = 1:size(SPM,1)
        SPM_M = [SPM_M;SPM{kk,rr}]; 
        S4M_M = [S4M_M;S4M{kk,rr}];
    end
    %unit conversion from cycles to radians for sigmaphi values
    if ~isempty(SPM_M)   
        SPM_M(:,2) = SPM_M(:,2)*2*pi;
        maxsp = max(SPM_M(:,2));
        maxs4 = max(S4M_M(:,2));
    else
        maxsp = 0;
        maxs4 = 0;
    end
    MAXSP = [MAXSP;maxsp];
    MAXS4 = [MAXS4;maxs4];
    SPM_elm = SPM_M;
    SP4_elm = S4M_M;
%     if strcmp(case_folder(end-4:end-1),'pfrr')
%         %folder_path for 2013 Poker Flat data 
%         op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%     else
%         %folder_path for 2013 Calgary data
%         op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);    
%     end       
    [~,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
    outfilename = strcat('prn_files_',signal,sep,'spdata.mat');
    save([op_path,outfilename],'SPM_elm','SP4_elm');
    %further eliminate lower sigmaphi values with sigmaphi mask
    if ~isempty(SPM_M)
        SPM_M = SPM_M(SPM_M(:,2)>spmask,:);
        S4M_M = S4M_M(S4M_M(:,2)>s4mask,:);
    else
        SPM_M = [];
        S4M_M = [];
    end
    outfilename = strcat('prn_files_',signal,sep,'fltrdata.mat');
    save([op_path,outfilename],'SPM_M','S4M_M');
end

MAXSP
MAXS4
splim = max(MAXSP);
s4lim = max(MAXS4);
fig_name = strcat(signal,'_LowRateSigmaPhiWithElevationMask_');
fig_1 = strcat(op_path,fig_name,'_fig1','.eps');
fig_2 = strcat(op_path,fig_name,'_fig2','.eps');
saveas(fig1,fig_1,'epsc2');
saveas(fig2,fig_2,'epsc2');

