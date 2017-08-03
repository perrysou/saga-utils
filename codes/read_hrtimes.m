function [hrtimes] = read_hrtimes(hrtimes_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
hrdata_id = fopen(hrtimes_path);
hrdata = textscan(hrdata_id,'%f %f %f %2f%2f %2f%2f') ;
prnlist = hrdata{1};
hrtimes = [];
for kk = 1:size(prnlist,1)
    prn = prnlist(kk,:);
    doy = hrdata{3}(kk);
    year = hrdata{2}(kk);
    tstthh = hrdata{4}(kk);
    tsttmm = hrdata{5}(kk);
    tendhh = hrdata{6}(kk);
    tendmm = hrdata{7}(kk);    
    dateV = datevec(datenum([year zeros(1,5)])+doy);
    tstt = datenum([dateV(1:3) tstthh tsttmm 0]);
    tend = datenum([dateV(1:3) tendhh tendmm 0]);  
    tspan = [tstt tend];
    hrtimes = [hrtimes;prn tspan];
end
end

