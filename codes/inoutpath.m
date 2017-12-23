function [in_path,op_path] = inoutpath(cases_folder,home_dir,year,doy,rcvr_name)
%specify the input and output path
%cases_foler: where local data is located
%home_dir: where processed data is located
%8.7 new data structure looks like this
%'/data1/public/Data/cases/pfrr/2013/342/grid108/'
%old data structure looks like this
%'/data1/public/cases/pfrr/grid108/2013/342/'
sep = filesep;

if strcmp(cases_folder(end-4:end-1),'pfrr')
    %folder_path for 2013 Poker Flat data 
    in_path = [cases_folder,year,sep,doy,sep,rcvr_name,sep]; 
    op_path = [home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep];
elseif strcmp(cases_folder(end-5:end-1),'tests')
    in_path = [cases_folder,year,sep,doy,sep,rcvr_name,sep]; 
    op_path = [home_dir,'mat/',rcvr_name,sep,year,sep,doy,sep];
else
    %folder_path for 2013 Calgary data
    year = 2013;
    in_path = [cases_folder,year,sep,doy,sep,rcvr_name,sep];
    op_path = [home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep];  
end

end

