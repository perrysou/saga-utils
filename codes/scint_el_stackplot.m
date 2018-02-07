function [MSP,MS4,FLTRDATA,CST,rcvr_op,tlim,splim,s4lim,signal] = scint_el_stackplot(cases_folder,home_dir,signal_type,doy,year,spmask,s4mask)
sep = filesep;
%specify signal type
switch signal_type
   case 0
      signal = 'L1CA';  %// GPS L1 legacy civil code 
   case 1
      signal = 'L2CM';  %// GPS L2 civil L code
   case 2
      signal = 'L2CL';  %// GPS L2 M+L combined tracking
   case 3
      signal = 'L2CLM';  %// GPS L5 civil in-phase
   case 4
      signal = 'L5I';   %// GPS L5 civil quadrature
   case 5
      signal = 'L5Q';
   case 6
      signal = 'L5IQ';
   case 7
      signal = 'L1CA-ALT1';
   case 8
      signal = 'CDMA-UHF-PILOT';
   case 9
      signal = 'CDMA-UHF-SYNC';
   otherwise
      error('Unknown signal.')
end
%specify path
% case_folder = '/media/Misc/MMAE552/~sdattaba/downloads/';
% case_folder = '/home/yang/cases/';
REF = [];
%operational receivers structure for the day
rcvr_op = [];
%receiver structure including all poker flat receivers
rcvr_struct = rx_dirs(cases_folder,year,doy);
if strcmp(year,'2015') && (strcmp(doy,'076') || strcmp(doy,'077'))
    rcvr_struct = ['ASTRArx';rcvr_struct;]
end
%% Read data from receiver structure
for rr = 1:size(rcvr_struct,1)  
    rcvr_name = rcvr_struct(rr,:)
%         if strcmp(cases_folder(end-4:end-1),'pfrr')
%             %folder_path for 2013 Poker Flat data 
%             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%             in_path = strcat([cases_folder,rcvr_name,sep,year,sep,doy,sep]);
%         else
%             %folder_path for 2013 Calgary data
%             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);  
%             in_path = strcat([cases_folder,rcvr_name,sep,doy,sep]);
%             year = 2013;
%         end   
    %in_path0: /data1/public/Data/cases/pfrr/
    
    offsett_str = datestr(datenum([str2num(year),1,0,0,0,0])+(str2num(doy))...
        -3600/24/3600,'mm/dd/yyyy HH:MM:SS')
    [~,offsett] = system(['date -ud ''',offsett_str,''' +%Y%j%H']);
    offsetyear = offsett(1:4);
    offsetdoy = offsett(5:7);
    offsethour = offsett(end-2:end-1);
%     [in_path_offset,~] = inoutpath(cases_folder,home_dir,offsetyear,offsetdoy,rcvr_name);
    
    [in_path_offset,~,in_path_tmp_offset] = sub_read_file(cases_folder,home_dir,offsetyear,offsetdoy,rcvr_name,sep)
    [in_path,op_path,in_path_tmp] = sub_read_file(cases_folder,home_dir,year,doy,rcvr_name,sep)
    if isempty(in_path) && isempty(in_path_tmp)
        continue;
    end
%     [in_path0,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
%     flag0 = dir([in_path0,'bin',sep,'*',year,'_',doy,'_*.bin']);
%     %in_path1: /data1/public/Data/cases/pfrr/from_usb/
%     in_path1 = [cases_folder,'from_usb',sep,year,sep,doy,sep,rcvr_name,sep]; 
%     flag1 = dir([in_path1,'bin',sep,'*',year,'_',doy,'_*.bin']);
%     %in_path2: /data1/from_usb/
%     in_path2 = ['/data1/from_usb',sep,year,sep,doy,sep,rcvr_name,sep];
%     flag2 = dir([in_path2,'bin',sep,'*',year,'_',doy,'_*.bin']);
%     if isempty(flag0) && isempty(flag1) && isempty(flag2)
%         disp(['Unable to find any binaries for ',rcvr_name,' on the server, please download manually']);
%         continue;
%     elseif ~isempty(flag0)
%         in_path = in_path0;
%     elseif ~isempty(flag1)
%         in_path = in_path1;
%     elseif ~isempty(flag2)
%         in_path = in_path2; 
%     end
%     in_path
%     in_path_tmp = in_path;
%     %repack and unpack processing
%     
% %     %repack first to remove any log files during uncompleted unpacking
% %     repack_comm = ['sh ',home_dir,'repack_v1.sh ',in_path,' ',op_path];
% %     system(repack_comm);
%     
%     %if binaries have been already unpacked (before automatic on-demand unpacking)
%     logflagout = dir([op_path,'txt',sep,'*',year,'_',doy,'_*.log']);
%     logflagin = dir([in_path,'txt',sep,'*',year,'_',doy,'_*.log']);
% %     keyboard;
%     if ~isempty(logflagin) && strcmp(in_path,in_path0)
%         disp(['binaries already unpacked;locate logfiles in ',in_path,'txt',sep]);
%     elseif isempty(logflagout)
%         %unpack 
%         unpack_comm = ['sh ',home_dir,'unpack_v1.sh ',in_path,' ',op_path]
%         system(unpack_comm); 
%         in_path = op_path;
%     else
%         disp(['binaries already unpacked;locate logfiles in ',op_path,'txt',sep]);
%         in_path = op_path;
%     end 
    rcvr_op = [rcvr_op; rcvr_name];  
    %Create output data folders
    command = strcat('mkdir -p ',{' '},op_path,'prn_files_',signal);
    system(cell2mat(command));
    %read_data
    [DATA, DATA_el,ref] = read_data(offsetyear,offsetdoy,in_path_offset,...
        doy,year,in_path,op_path,sep,signal_type,signal);
%         size(DATA)
    REF = [REF;ref];
    read_navsol(offsetyear,offsetdoy,in_path_offset,...
        doy,year,sep,in_path,op_path,signal);
    read_iono(doy,year,sep,in_path,op_path,signal);
    read_channel(doy,year,sep,in_path,op_path,signal);
    %repacking
    %2017/9/11 commented out repack commands just to be safe
%     repack_comm = ['sh ',home_dir,'repack_v1.sh ',in_path_tmp,' ',op_path];
%     system(repack_comm);
%     repack_comm = ['sh ',home_dir,'repack_v1.sh ',in_path_tmp_offset,' ',op_path];
%     system(repack_comm);
end
refprn = unique(REF);
if isempty(rcvr_op)
    disp(['Unable to continue w/o any logfile of any receiver for doy ',...
        num2str(doy,'%03i'),'.  Now exiting...']);
    MSP = [];FLTRDATA = [];CST = [];tlim = [];splim = [];MS4 = [];s4lim = [];
    return;
else
%% Plot elevation masked(el>30 deg) sigmaphi and save masked sigmaphi(with el>30 deg)
[tlim,splim,s4lim] = el_and_masked_sigmaphi(cases_folder,home_dir,rcvr_op,doy,year,signal,spmask,s4mask);
close all;

%% Plot masked sigmaphi(with el>30 deg) before finding common scintillation times
MSP = [];
MS4 = [];
for kk = 1:32
    if kk<=16
    fig3 = subplot(8,2,kk);   
    elseif kk == 17
    figure;
    subplot(8,2,kk-16);
    else
    fig4 = subplot(8,2,kk-16);
    end
    RCVRNAME = [];
    for rr = 1:size(rcvr_op,1)
        rcvr_name = rcvr_op(rr,:);
%         if strcmp(cases_folder(end-4:end-1),'pfrr')
%             %folder_path for 2013 Poker Flat data 
%             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%         else
%             %folder_path for 2013 Calgary data
%             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);    
%         end
        [~,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
        outfilename = strcat('prn_files_',signal,sep,'fltrdata.mat');
        load([op_path,outfilename]);
        if ~isempty(SPM_M)
            svdata = SPM_M(SPM_M(:,3)==kk,:);
            if ~isempty(svdata)
                time = svdata(:,1);
                sigmaphi = svdata(:,2);
                MSP = [MSP;time sigmaphi ones(size(time))*kk ones(size(time))*rr];
                dt = 1/24/3600*600;
                [time_e,sigmaphi_e] = discont_proc(time,sigmaphi,dt);
                [color] = rx_color(rcvr_name);
                plot(time_e,sigmaphi_e,'LineWidth',1,'Color',color);  
                RCVRNAME = [RCVRNAME;rcvr_name];
                hold on;
            end
        end
        if ~isempty(S4M_M)
            svdata_s4 = S4M_M(S4M_M(:,3)==kk,:);
            if ~isempty(svdata_s4)
                time_s4 = svdata_s4(:,1);
                s4 = svdata_s4(:,2);
                MS4 = [MS4;time_s4 s4 ones(size(time_s4))*kk ones(size(time_s4))*rr];
                dt = 1/24/3600*600;
                [time_e,s4_e] = discont_proc(time_s4,s4,dt);
%                 [color] = rx_color(rcvr_name);
%                 plot(time_e,s4_e,'LineWidth',1,'Color',color);  
%                 RCVRNAME = [RCVRNAME;rcvr_name];
%                 hold on;
            end
        end
    end
    set(gca,'FontSize',8,'YTick',[0 splim],...,
            'YTickLabel',{'0','MAX'},'YGrid','on');
    if isempty(tlim)
        datetick('x','HH','keeplimits');
    else
        xlim(tlim);
        datetick('x','HH','keeplimits');
    end
    ylim([0 splim*1.2]);
    xlabel(['PRN',num2str(kk)]); 
    if kk == 32 || kk == 16
        ylabel('$\sigma_\phi$ in [rad]');
        legend(RCVRNAME);
    end
end
fig_name = strcat(signal,'_MaskedSigmaphi');
fig_3 = strcat(op_path,fig_name,'_fig3','.eps');
fig_4 = strcat(op_path,fig_name,'_fig4','.eps'); 
saveas(fig3,fig_3,'epsc2');
saveas(fig4,fig_4,'epsc2');
close all;

%% Plot sigmaphi over common scintillation times
FLTRDATA = [];
CST = [];
for kk = 1:32
    if kk<=16
    fig5 = subplot(8,2,kk);   
    elseif kk == 17
    figure;
    subplot(8,2,kk-16);
    else
    fig6 = subplot(8,2,kk-16);
    end
    RCVRNAME = [];   
    for rr = 1:size(rcvr_op,1)        
        rcvr_name = rcvr_op(rr,:);
%         if strcmp(cases_folder(end-4:end-1),'pfrr')
%             %folder_path for 2013 Poker Flat data 
%             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%         else
%             %folder_path for 2013 Calgary data
%             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);    
%         end
        [~,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
        load([op_path,strcat('prn_files_',signal,sep,'fltrdata.mat')]);
        if ~isempty(SPM_M)
        svdata = SPM_M(SPM_M(:,3)==kk,:);
        if ~isempty(svdata)
            time = svdata(:,1);
            tdiff = diff(time);
            disc = find(tdiff>dt);
            if size(time,1)<=2
            trx{rr} = sortrows([time(1);time(end)]);
            else
            trx{rr} = sortrows([time(1);time(disc);time(disc+1);time(end)]);
            end
            trx{rr} = trx{rr}';
        else trx{rr} = [];
        end
        else trx{rr} = [];
        end
    end
    %find general scintillation times
    t = find_common_times(trx); 
    if ~isempty(t)
        CST = [CST; t ones(size(t))*kk];
        for rr = 1:size(rcvr_op,1)
            SIGMAPHI = [];TIME = [];          
            rcvr_name = rcvr_op(rr,:);
%             if strcmp(cases_folder(end-4:end-1),'pfrr')
%                 %folder_path for 2013 Poker Flat data 
%                 op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
%             else
%                 %folder_path for 2013 Calgary data
%                 op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);    
%             end
            [~,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
            load([op_path,strcat('prn_files_',signal,sep,'spdata.mat')]);
            if ~isempty(SPM_elm)
            svdata = SPM_elm(SPM_elm(:,3)==kk,:);
            time = svdata(:,1);
            if ~isempty(svdata)
                for tt = 1:2:size(t,1)
                    svdata_masked = svdata(time<=t(tt+1) & time>=t(tt),:);
                    time_d = svdata_masked(:,1);
                    sigmaphi_d = svdata_masked(:,2);
                    SIGMAPHI = [SIGMAPHI;sigmaphi_d];
                    TIME = [TIME;time_d];
                end
                fltrdata = [TIME SIGMAPHI ones(size(TIME))*kk ones(size(TIME))*rr];
                FLTRDATA = [FLTRDATA;fltrdata];
                [time_e,sigmaphi_e] = discont_proc(TIME,SIGMAPHI,dt); 
                [color] = rx_color(rcvr_name);
                plot(time_e,sigmaphi_e,'LineWidth',1,'Color',color);
                RCVRNAME = [RCVRNAME;rcvr_name];
                hold on;  
            end
            end
        end
    end
    
    set(gca,'FontSize',8,'YTick',[0 splim],...,
            'YTickLabel',{'0','MAX'},'YGrid','on');
    if isempty(tlim)
        datetick('x','HH','keeplimits');
    else
        xlim(tlim);
        datetick('x','HH','keeplimits');
    end
    ylim([0 splim*1.2]);
    xlabel(['PRN',num2str(kk)]);
    
    if kk == 32 || kk == 16
        ylabel('$\sigma_\phi$ in [rad]');
        legend(RCVRNAME);
    end
end
    
fig_name = strcat(signal,'_CommonScintillationTimes');
fig_5 = strcat(op_path,fig_name,'_fig5','.eps');
fig_6 = strcat(op_path,fig_name,'_fig6','.eps');
saveas(fig5,fig_5,'epsc2');
saveas(fig6,fig_6,'epsc2');
close all;
% MSP = sortrows(MSP,[3 1]);
% % MISSDATA = setdiff(MSP,FLTRDATA,'rows');
% % MISSDATA = sortrows(MISSDATA,[3 1]);
% % HITDATA = sortrows(FLTRDATA,[3 1]);
end
end

function [in_path,op_path,in_path_tmp] = sub_read_file(cases_folder,home_dir,year,doy,rcvr_name,sep)
[in_path0,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name);
flag0 = dir([in_path0,'bin',sep,'*',year,'_',doy,'_*.bin']);
%in_path1: /data1/public/Data/cases/pfrr/from_usb/
in_path1 = [cases_folder,'from_usb',sep,year,sep,doy,sep,rcvr_name,sep];
flag1 = dir([in_path1,'bin',sep,'*',year,'_',doy,'_*.bin']);
%in_path2: /data1/from_usb/
in_path2 = ['/data1/from_usb',sep,year,sep,doy,sep,rcvr_name,sep];
if strcmp(rcvr_name,'ASTRArx')
    flag2 = dir([in_path2,'bin',sep,'*.bin']);
else
    flag2 = dir([in_path2,'bin',sep,'*',year,'_',doy,'_*.bin']);
end
if isempty(flag0) && isempty(flag1) && isempty(flag2)
    disp(['Unable to find any binaries for ',rcvr_name,' on the server, please download manually']);
    in_path = [];in_path_tmp = [];
    return;
elseif ~isempty(flag0) && any([flag0.bytes])
    in_path = in_path0;
elseif ~isempty(flag1) && any([flag1.bytes])
    in_path = in_path1;
elseif ~isempty(flag2) && any([flag2.bytes])
    in_path = in_path2; 
else
    disp(['Unable to find any binaries for ',rcvr_name,' on the server, please download manually']);
    in_path = [];in_path_tmp = [];
    return;
end
in_path;
in_path_tmp = in_path;
%repack and unpack processing

%     %repack first to remove any log files during uncompleted unpacking
%     repack_comm = ['sh ',home_dir,'repack_v1.sh ',in_path,' ',op_path];
%     system(repack_comm);

%if binaries have been already unpacked (before automatic on-demand unpacking)
if strcmp(rcvr_name,'ASTRArx')
    logflagout = dir([op_path,'txt',sep,'*.log']);
    logflagin = dir([in_path,'txt',sep,'*.log']);
else
    logflagout = dir([op_path,'txt',sep,'*',year,'_',doy,'_*.log']);
    logflagin = dir([in_path,'txt',sep,'*',year,'_',doy,'_*.log']);
end

if ~isempty(logflagin) && strcmp(in_path,in_path0)
    disp(['binaries already unpacked;locate logfiles in ',in_path,'txt',sep]);
elseif isempty(logflagout)
    %unpack 
    unpack_comm = ['sh ',home_dir,'unpack_v1.sh ',in_path,' ',op_path]
    system(unpack_comm); 
    in_path = op_path;
else
    disp(['binaries already unpacked;locate logfiles in ',op_path,'txt',sep]);
    in_path = op_path;
end 
end