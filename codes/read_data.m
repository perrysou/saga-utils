function [DATA,DATA_el,REFPRN] = read_data(offsetyear,offsetdoy,in_path_offset,...
    doy,year,folder_path,op_path,sep,signal_type,signal)
%% scint
%Read and combine scint files
if strcmp(folder_path(end-7:end-1),'ASTRArx') || strcmp(folder_path(end-16:end-10),'ASTRArx')
    scintfilestruct = dir([folder_path,'txt',sep,'scint*.log']);
    scintfilestruct0 = [];
else
    scintfilestruct = dir([folder_path,'txt',sep,'scint*',year,'_',doy,'_*.log']);
    scintfilestruct0 = dir([in_path_offset,'txt',sep,'scint*',offsetyear,'_',offsetdoy,'_*.log']);
end

if ~isempty(scintfilestruct0)
    scintfilename0 = dlmread([in_path_offset,'txt',sep,scintfilestruct0(end).name]);
    [in_path_offset,'txt',sep,scintfilestruct0(end).name]
else 
    scintfilename0 = [];
end

if ~isempty(scintfilestruct)
    SCINT=[];
    for jj=1:size(scintfilestruct,1)
        scintfile = scintfilestruct(jj,:);
        scintfilename = dlmread(strcat([folder_path,'txt',sep,scintfile.name]));
        st = datestr(gps2utc(scintfilename(1,1:2)),'HHMM-');
        se = datestr(gps2utc(scintfilename(end,1:2)),'HHMM UT');
        disp([scintfile.name,' actually has data for ',st,se]);
        SCINT = [SCINT;scintfilename];
    end
    SCINT = [SCINT;scintfilename0];
    
    size(find(SCINT(:,1)>3640));
    SCINT = SCINT(SCINT(:,1)<=3640,:);
    size(find(SCINT(:,1)>3640));
    SCINT = sortrows(SCINT,[1 2 3]);
    %output scintdata for the whole day
    scfile = SCINT;
    
    %Specify signal type
    SCINT = SCINT((SCINT(:,14)==signal_type),:);
    %Remove invalid entries
    SCINT = SCINT((SCINT(:,8)<1)&SCINT(:,8)>-1&SCINT(:,5)>-1,:);
    REFPRN = unique(SCINT((SCINT(:,13)==1|SCINT(:,13)==2),15));
    %Gather data
    ORTW = SCINT(:,1);
    ORTS = SCINT(:,2)+SCINT(:,3);
    S4 = SCINT(:,5);
    SIGMAPHI = SCINT(:,8);
    TAU0 = SCINT(:,11);
    
    % SIGMAPHI = SCINT(:,8);
    PRN = SCINT(:,15);

    DATA = [ORTW ORTS SIGMAPHI PRN S4 TAU0];
%     gps2utc(DATA(1:10   ,1:2))
%     keyboard;
    % DATA = sortrows(DATA,4);
else
    DATA = [];
    REFPRN = [];
    scfile = [];
end
outfilename = strcat('prn_files_',signal,sep,'scintdata.mat');
command = strcat('mkdir -p ', [' ', op_path, 'prn_files_', signal, sep]);
system(command);
save([op_path,outfilename],'DATA','scfile');
% keyboard;
%% txinfo
%Read and combine txinfo files
if strcmp(folder_path(end-7:end-1),'ASTRArx') || strcmp(folder_path(end-16:end-10),'ASTRArx')
    azelfilestruct = dir([folder_path,'txt',sep,'txinfo*.log']);
    azelfilestruct0 = [];
else
    azelfilestruct = dir([folder_path,'txt',sep,'txinfo*',year,'_',doy,'_*.log']);
    azelfilestruct0 = dir([in_path_offset,'txt',sep,'txinfo*',offsetyear,'_',offsetdoy,'_*.log']);
end

if ~isempty(azelfilestruct0)
    azelfilename0 = dlmread([in_path_offset,'txt',sep,azelfilestruct0(end).name]);
    [in_path_offset,'txt',sep,azelfilestruct0(end).name]
else 
    azelfilename0 = [];
end

if ~isempty(azelfilestruct)
    AZEL=[];
    for jj=1:size(azelfilestruct,1)
        azelfile = azelfilestruct(jj,:);
        azelfilename = dlmread(strcat([folder_path,'txt',sep,azelfile.name]));
        st = datestr(gps2utc(azelfilename(1,1:2)),'HHMM-');
        se = datestr(gps2utc(azelfilename(end,1:2)),'HHMM UT');
%         disp([azelfile.name,' actually has data for ',st,se]);
        AZEL = [AZEL;azelfilename];
    end
    AZEL = [AZEL;azelfilename0];
    
    size(find(AZEL(:,1)>3640))
    AZEL = AZEL(AZEL(:,1)<=3640,:);
    size(find(AZEL(:,1)>3640))
    AZEL = sortrows(AZEL,[1 2 3]);
    %Remove invalid entries
    AZEL = AZEL((AZEL(:,5)>0),:);
    AZEL_old = AZEL;
    %Low elevation threshold;
    AZEL_low = AZEL(AZEL(:,5)>30,:);
    % AZEL_low = sortrows(AZEL_low,8);

    %Gather data
    ORTW_EL = AZEL(:,1);
    ORTS_EL = AZEL(:,2)+AZEL(:,3);
    AZ = AZEL(:,4);
    EL = AZEL(:,5); 
    PRN_EL = AZEL(:,8);

    DATA_el = [ORTW_EL ORTS_EL EL PRN_EL AZ];
    % DATA = sortrows(DATA,4);
else
    DATA_el = [];
end
outfilename = strcat('prn_files_',signal,sep,'txinfodata.mat');
command = strcat('mkdir -p ', [' ', op_path, 'prn_files_', signal, sep]);
system(command);
save([op_path,outfilename],'DATA_el');

end

